#include "err.h"
#include "fp.h"
#include "sigs.h"
#include "size.h"
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Static Constants (N) -------------------------------------------------------

#define NDIM ((size_t) 1)
#define NPOLY ((size_t) 1)
#define NBEAD ((size_t) 64)
#define NTSTEP ((size_t) 1 << 14)
#define NPSTEP ((size_t) 1 << 18)
#define NTRSTEP ((size_t) 1 << 4)
#define NPRSTEP ((size_t) 1 << 8)
#define NSUBDIV ((size_t) 16)

/*
These assertions guarantee that adding or multiplying two steps or indices
never wraps around, which makes their manipulation easier.
*/
static_assert(NTSTEP <= SQRT_SIZE_MAX, "too many thermalization steps");
static_assert(NPSTEP <= SQRT_SIZE_MAX, "too many production steps");
static_assert(NTRSTEP <= NTSTEP, "too many thermalization recording steps");
static_assert(NPRSTEP <= NPSTEP, "too many production recording steps");

// Types () -------------------------------------------------------------------

struct bead {
  double d[NDIM];
};

struct poly {
  size_t from;
  size_t to;
  struct bead r[NBEAD];
};

struct napkin {
  gsl_rng* rng; // State of the random number generator.
  struct bead origin;
  struct {
    struct {
      size_t accepted;
      size_t rejected;
      double dx;
    } ss; // Single slice.
    struct {
      size_t accepted;
      size_t rejected;
      double dx;
    } comd; // Center of mass displacement.
  } params; // Optimization parameters of each move.
  union {
    struct {
      size_t ipoly;
      size_t ibead;
      struct bead r;
    } ss;
    struct {
      size_t ipoly;
      struct poly R;
    } comd;
  } history; // Undo history of the latest move.
  void (* accept)(void); // Procedures for processing the latest move.
  void (* reject)(void);
  void (* adjust)(void);
  size_t itstep; // Thermalization steps taken.
  size_t ipstep; // Production steps taken.
  size_t itrstep; // Recording steps taken.
  size_t iprstep; // Recording steps taken.
  double L; // Lattice constant.
  double beta;
  double lambda;
  double tau;
  double (* Vext)(struct bead); // External potential.
  double (* V)(struct bead); // Potential energy.
  double (* K)(struct bead); // Kinetic energy.
  double (* Vend)(struct bead); // Potential energy at open ends.
  double (* Kend)(struct bead); // Kinetic energy at open ends.
  struct {
    size_t N; // Samples.
    double mean; // Mean.
    double M2; // Second moment.
    double var; // Variance.
    double std; // Standard deviation.
    double stderr; // Standard error of the mean.
  } tde; // Thermodynamic estimator.
  struct poly R[NPOLY];
};

// Lost Souls () --------------------------------------------------------------

static size_t ran_index(gsl_rng* const rng, size_t const n) {
  return (size_t) gsl_rng_uniform_int(rng, n);
}

// Global State () ------------------------------------------------------------

static struct napkin napkin;

// Vector Math (bead) ---------------------------------------------------------

static struct bead bead_uwrap(struct bead const r) {
  struct bead s;

  for (size_t idim = 0; idim < NDIM; ++idim)
    s.d[idim] = fp_uwrap(r.d[idim], napkin.L);

  return s;
}

static struct bead bead_wrap(struct bead const r) {
  struct bead s;

  for (size_t idim = 0; idim < NDIM; ++idim)
    s.d[idim] = fp_wrap(r.d[idim], napkin.L);

  return s;
}

// This only takes the closest periodic image into account.
static struct bead bead_sub(struct bead const r0, struct bead const r1) {
  struct bead r;

  for (size_t idim = 0; idim < NDIM; ++idim)
    r.d[idim] = r1.d[idim] - r0.d[idim];

  return bead_wrap(r);
}

static double bead_norm2(struct bead const r) {
  double s = 0;

  for (size_t idim = 0; idim < NDIM; ++idim)
    s += gsl_pow_2(r.d[idim]);

  return s;
}

static double bead_norm(struct bead const r) {
  return sqrt(bead_norm2(r));
}

// Configuration Forcing (force) ----------------------------------------------

/*
Close the current configuration.
*/
static void force_close(void) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    napkin.R[ipoly].from = ipoly;
    napkin.R[ipoly].to = ipoly;
  }
}

/*
Open the current configuration.
*/
static void force_open(void) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    napkin.R[ipoly].from = SIZE_MAX;
    napkin.R[ipoly].to = SIZE_MAX;
  }
}

// Initial Configurations (conf) ----------------------------------------------

/*
Random initial configuration.
*/
static void conf_rand(void) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      for (size_t idim = 0; idim < NDIM; ++idim)
        napkin.R[ipoly].r[ibead].d[idim] =
          gsl_ran_flat(napkin.rng, 0, napkin.L);
}

/*
Random point initial configuration.
*/
static void conf_randpt(void) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t idim = 0; idim < NDIM; ++idim) {
      double const x = gsl_ran_flat(napkin.rng, 0, napkin.L);

      for (size_t ibead = 0; ibead < NBEAD; ++ibead)
        napkin.R[ipoly].r[ibead].d[idim] = x;
    }
}

/*
Random lattice initial configuration.
*/
static void conf_randlatt(void) {
  double const w = ceil(pow(NPOLY, 1.0 / NDIM));
  double const v = napkin.L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
      size_t m = ipoly;

      for (size_t idim = 0; idim < NDIM; ++idim) {
        size_div_t const z = size_div(m, n);

        napkin.R[ipoly].r[ibead].d[idim] =
          (double) z.rem * v + gsl_ran_flat(napkin.rng, 0, v);
        m = z.quot;
      }
    }
}

/*
Point lattice initial configuration.
*/
static void conf_ptlatt(void) {
  double const w = ceil(pow(NPOLY, 1.0 / NDIM));
  double const v = napkin.L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
      size_t m = ipoly;

      for (size_t idim = 0; idim < NDIM; ++idim) {
        size_div_t const z = size_div(m, n);

        napkin.R[ipoly].r[ibead].d[idim] = (double) z.rem * v;
        m = z.quot;
      }
    }
}

/*
Circle lattice initial configuration.
The equations are based on the theory of Lissajous knots.
*/
static void conf_circlatt(void) {
  double const w = ceil(pow(NPOLY, 1.0 / NDIM));
  double const v = napkin.L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
      double const t = (double) ibead / NBEAD;

      size_t m = ipoly;

      for (size_t idim = 0; idim < NDIM; ++idim) {
        double const phi = (double) idim / NDIM;
        double const y = sin(M_2PI * (t + phi / 2)) / 2;

        size_div_t const z = size_div(m, n);
        napkin.R[ipoly].r[ibead].d[idim] = ((double) z.rem + (y + 1) / 2) * v;
        m = z.quot;
      }
    }
}

// Physical Moves (move) ------------------------------------------------------

/*
r <=> a: z >=< 1
r = 0: z = 3 / 2
a = 0: z = 1 / 2
*/
static double surface(double const a, double const r) {
  return 0.5 + a / (a + r);
}

/*
r <=> a: z >=< 1
r = 0: z = 1 + exp(-a)
a = 0: z = 1 - exp(-r)
*/
static double another_surface(double const a, double const r) {
  return 1 - exp(-a) + exp(-r);
}

static void move_accept_ss(void) {
  ++napkin.params.ss.accepted;
}

static void move_reject_ss(void) {
  napkin.R[napkin.history.ss.ipoly].r[napkin.history.ss.ibead] =
    napkin.history.ss.r;

  ++napkin.params.ss.rejected;
}

static void move_adjust_ss(void) {
  napkin.params.ss.dx = napkin.params.ss.dx *
    surface((double) napkin.params.ss.accepted,
        (double) napkin.params.ss.rejected);
}

static void move_ss(size_t const ipoly, size_t const ibead) {
  napkin.accept = move_accept_ss;
  napkin.reject = move_reject_ss;
  napkin.adjust = move_adjust_ss;

  napkin.history.ss.ipoly = ipoly;
  napkin.history.ss.ibead = ibead;
  napkin.history.ss.r = napkin.R[ipoly].r[ibead];

  for (size_t idim = 0; idim < NDIM; ++idim)
    napkin.R[ipoly].r[ibead].d[idim] = fp_uwrap(
        napkin.R[ipoly].r[ibead].d[idim] +
        napkin.params.ss.dx * gsl_ran_flat(napkin.rng, -1, 1),
        napkin.L);
}

static void move_accept_comd(void) {
  ++napkin.params.comd.accepted;
}

static void move_reject_comd(void) {
  napkin.R[napkin.history.comd.ipoly] = napkin.history.comd.R;

  ++napkin.params.comd.rejected;
}

static void move_adjust_comd(void) {
  napkin.params.comd.dx = napkin.params.comd.dx *
    another_surface((double) napkin.params.comd.accepted,
        (double) napkin.params.comd.rejected);
}

static void move_comd(size_t const ipoly) {
  napkin.accept = move_accept_comd;
  napkin.reject = move_reject_comd;
  napkin.adjust = move_adjust_comd;

  napkin.history.comd.ipoly = ipoly;
  napkin.history.comd.R = napkin.R[ipoly];

  for (size_t idim = 0; idim < NDIM; ++idim) {
    double const x = napkin.params.comd.dx * gsl_ran_flat(napkin.rng, -1, 1);

    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      napkin.R[ipoly].r[ibead].d[idim] = fp_uwrap(
          napkin.R[ipoly].r[ibead].d[idim] + x,
          napkin.L);
  }
}

static void move_accept_bisect(void) {
  err_abort(NULL);
}

static void move_reject_bisect(void) {
  err_abort(NULL);
}

static void move_bisect(size_t const ipoly, size_t const ibead) {
  err_abort(NULL);
}

// Constants () ---------------------------------------------------------------

static double const hbar = 1;
static double const m = 1;
static double const k = 1;
static double const epsilon = 1; // L--J 6--12 P
static double const sigma = 1; // L--J 6--12 P
static double const omega = 1; // HP

// System-Specific Quantities () ----------------------------------------------

static double Kbead(size_t const ipoly, size_t const ibead) {
  double K = 0;

  for (int i = 0; i < 2; ++i) {
    bool const p = ibead == 0 && i == 0;
    bool const q = ibead == NBEAD - 1 && i == 1;
    double (* const f)(struct bead) = p || q ? napkin.Kend : napkin.K;
    size_t const jbead = size_uwrap(i == 0 ? ibead - 1 : ibead + 1, NBEAD);
    size_t const jpoly = p ? napkin.R[ipoly].from :
      q ? napkin.R[ipoly].to :
      ipoly;

    if (jpoly != SIZE_MAX)
      K += f(bead_sub(napkin.R[ipoly].r[ibead], napkin.R[jpoly].r[jbead]));
  }

  return K;
}

static double Kpoly(size_t const ipoly) {
  double K = 0;

  for (size_t ibead = 0; ibead < NBEAD; ++ibead)
    K += Kbead(ipoly, ibead);

  return K / 2;
}

// The log(2 * M_2PI * lambda * tau) term is not here
// due to scaling and offsets?
static double Ktotal(void) {
  double K = 0;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    K += Kpoly(ipoly);

  return K;
}

static double Vextbead(size_t const ipoly, size_t const ibead) {
  return napkin.Vext(napkin.R[ipoly].r[ibead]);
}

static double Vbeads(size_t const ibead) {
  double V = 0;

  bool const p = ibead == 0 || ibead == NBEAD - 1;
  double (* const f)(struct bead) = p ? napkin.Vend : napkin.V;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < NPOLY; ++jpoly)
      V += f(bead_sub(napkin.R[ipoly].r[ibead], napkin.R[jpoly].r[ibead]));

  return V;
}

static double Vtotal(void) {
  double V = 0;

  for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
    V += Vbeads(ibead);

    for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
      V += Vextbead(ipoly, ibead);
  }

  return V;
}

static double Stotal(void) {
  return Vtotal() + Ktotal();
}

// Estimators (est) -----------------------------------------------------------

// \langle E_T\rangle & = \frac 1{M \tau} \sum_{k = 1}^M
// \langle \frac{d N} 2 - \frac{|R_{k + 1 \bmod M} - R_k|^2}{4 \lambda \tau} + \tau V(R_k)\rangle
static double est_tde(void) {
  double const KC = (double) NDIM / 2;
  double const K = Ktotal() / (4 * napkin.lambda * napkin.tau) / (NBEAD * NPOLY);
  double const V = Vtotal() * napkin.tau / (NBEAD * NPOLY);
  return (KC - K + V) / napkin.tau;
}

static void est_gather_tde(void) {
   double const E = est_tde();
   ++napkin.tde.N;
   double const Delta = E - napkin.tde.mean;
   napkin.tde.mean += Delta / (double) napkin.tde.N;
   napkin.tde.M2 += Delta * (E - napkin.tde.mean);

   napkin.tde.var = napkin.tde.M2 / (double) (napkin.tde.N - 1);
   napkin.tde.std = sqrt(napkin.tde.var);
   napkin.tde.stderr = sqrt(napkin.tde.var / (double) napkin.tde.N);
}

// Bad.
struct {
  size_t N; // Samples.
  double mean; // Mean.
  double M2; // Second moment.
  double var; // Variance.
  double std; // Standard deviation.
  double stderr; // Standard error of the mean.
} urgh; // Worldline histogram.
static void est_gather_x(void) {
   double const E = fp_wrap(napkin.R[0].r[NBEAD / 2].d[0], napkin.L);
   ++urgh.N;
   double const Delta = E - urgh.mean;
   urgh.mean += Delta / (double) urgh.N;
   urgh.M2 += Delta * (E - urgh.mean);
   urgh.var = urgh.M2 / (double) (urgh.N - 1);
   urgh.std = sqrt(urgh.var);
   urgh.stderr = sqrt(urgh.var / (double) urgh.N);
}

// Workers (work) -------------------------------------------------------------

static double work_ss(void) {
  size_t const ipoly = ran_index(napkin.rng, NPOLY);
  size_t const ibead = ran_index(napkin.rng, NBEAD);

  // Local `Vbeads`, because the rest stays the same.
  double const V0 = Vbeads(ibead) + Vextbead(ipoly, ibead);
  double const K0 = Kbead(ipoly, ibead);

  move_ss(ipoly, ibead);

  double const V1 = Vbeads(ibead) + Vextbead(ipoly, ibead);
  double const K1 = Kbead(ipoly, ibead);

  double const DeltaS_V = (V1 - V0) * napkin.tau;
  double const DeltaS_K = (K1 - K0) / (4 * napkin.lambda * napkin.tau);

  return DeltaS_V + DeltaS_K;
}

static double work_comd(void) {
  size_t const ipoly = ran_index(napkin.rng, NPOLY);

  // Global `Vtotal`, because `K` is unchanged (maybe not for permutations).
  double const V0 = Vtotal();

  move_comd(ipoly);

  double const V1 = Vtotal();

  double const DeltaS_V = (V1 - V0) * napkin.tau;

  return DeltaS_V;
}

// Work Horses () -------------------------------------------------------------

static void choose(double const DeltaS) {
  if (DeltaS <= 0)
    napkin.accept();
  else {
    double const q = gsl_ran_flat(napkin.rng, 0, 1);
    if (q <= exp(-DeltaS))
      napkin.accept();
    else
      napkin.reject();
  }

  napkin.adjust();
}

// Printing into Data Files (disp) --------------------------------------------

/*
Print bead into stream `fp`.
*/
static void disp_bead(FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  fprintf(fp, "%zu", iindex);

  for (size_t idim = 0; idim < NDIM; ++idim)
    fprintf(fp, " %f", napkin.R[ipoly].r[ibead].d[idim]);

  fprintf(fp, "\n");
}

/*
Print polymer into stream `fp`.
*/
static void disp_poly(FILE* const fp) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      disp_bead(fp, ibead, ipoly, ibead);

    if (napkin.R[ipoly].to != SIZE_MAX)
      disp_bead(fp, NBEAD, napkin.R[ipoly].to, 0);

    fprintf(fp, "\n");
  }
}

static void disp_tde(FILE* const fp) {
  fprintf(fp, "%zu %f %f %f\n",
      napkin.ipstep, est_tde(), napkin.tde.mean, napkin.tde.stderr * NBEAD);
}

/*
Print parameters into stream `fp`.
*/
static void disp_drift(FILE* const fp) {
  fprintf(fp, "%zu %f %f\n",
      napkin.itstep + napkin.ipstep,
      napkin.params.ss.dx, napkin.params.comd.dx);
}

// Saving into Data Files (save)

static void save_length(void) {
  FILE* const fp = fopen("qho-length.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  fprintf(fp, "%f\n", napkin.L);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

static void save_esl(void) {
  FILE* const fp = fopen("qho-ensemble.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  disp_poly(fp);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

// Other Crap () --------------------------------------------------------------

static void status_line(void) {
  printf("i / N = (%zu + %zu) / (%zu + %zu) = %zu %%\n"
      "E = %f +- %f\n"
      "K + V = %f + %f = %f\n",
      napkin.itstep, napkin.ipstep, NTSTEP, NPSTEP,
      100 * (napkin.itstep + napkin.ipstep) / (NTSTEP + NPSTEP),
      napkin.tde.mean, napkin.tde.stderr * NBEAD,
      Ktotal(), Vtotal(), Ktotal() + Vtotal());
}

static void status_file(void) {
  save_esl();
}

static double lj612b(struct bead const r) {
  double const sigmar2 = gsl_pow_2(sigma) / bead_norm2(r);

  return 4 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

static double harmb(struct bead const r) {
  return m * gsl_pow_2(omega) * bead_norm2(bead_wrap(r)) / 2;
}

static double sproink(void) {
  FILE* const fp = fopen("qho-subdivisions.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  fprintf(fp, "%zu\n", NSUBDIV);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

static double splat(void) {
  FILE* const fp = fopen("qho-potential.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  double const v = napkin.L / NSUBDIV;

  for (size_t i = 0; i < size_pow(NSUBDIV + 1, NDIM); ++i) {
    size_t m = i;

    struct bead r;

    for (size_t idim = 0; idim < NDIM; ++idim) {
      size_div_t const z = size_div(m, NSUBDIV + 1);

      r.d[idim] = (double) z.rem * v;
      m = z.quot;
    }

    for (size_t idim = 0; idim < NDIM; ++idim)
      fprintf(fp, "%f ", r.d[idim]);

    fprintf(fp, "%f\n", napkin.Vext(r));
  }

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

static double (* const workers[])(void) = {work_ss, work_comd};
#define NWORKERS (sizeof workers / sizeof *workers)

static void not_main(void) {
  err_reset();

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

  napkin.params.ss.accepted = 0;
  napkin.params.ss.rejected = 0;
  napkin.params.ss.dx = 1;
  napkin.params.comd.accepted = 0;
  napkin.params.comd.rejected = 0;
  napkin.params.comd.dx = 1;

  napkin.L = 10;
  double const T = 0.1;
  napkin.beta = 1 / (k * T);
  napkin.lambda = gsl_pow_2(hbar) / (2 * m);
  napkin.tau = napkin.beta / NBEAD;

  for (size_t idim = 0; idim < NDIM; ++idim)
    napkin.origin.d[idim] = 0;

  double const q = hbar * omega / 2;
  double const E = q / tanh(q * napkin.beta);
  printf("expect: E = %f\n", E);

  napkin.Vext = harmb;
  napkin.V = lj612b;
  napkin.Vend = lj612b;
  napkin.K = bead_norm2;
  napkin.Kend = bead_norm2;

  napkin.tde.N = 0;
  napkin.tde.mean = 0;
  napkin.tde.M2 = 0;
  napkin.tde.var = 0; // This should go elsewhere.
  napkin.tde.std = 0; // This should go elsewhere.
  napkin.tde.stderr = 0; // This should go elsewhere.

  gsl_rng_env_setup();
  napkin.rng = gsl_rng_alloc(gsl_rng_default);

  conf_circlatt();
  force_close();

  save_length();

  sproink();
  splat();

  FILE* const driftfp = fopen("qho-drift.data", "w");
  if (driftfp == NULL)
    err_abort(fopen);

  FILE* const tdefp = fopen("qho-tde.data", "w");
  if (tdefp == NULL)
    err_abort(fopen);

  // Just to prevent empty data files...
  disp_drift(driftfp);
  fflush(driftfp);
  fprintf(tdefp, "%zu %f %f %f\n", 0, E, E, 0);
  fflush(tdefp);

  napkin.itstep = 0;
  napkin.ipstep = 0;
  napkin.itrstep = 0;
  napkin.iprstep = 0;
  for (size_t istep = 0; istep < NTSTEP + NPSTEP; ++istep) {

    choose(workers[ran_index(napkin.rng, NWORKERS)]());

    int signum;
    if (sigs_use(&signum))
      switch (signum) {
        case SIGUSR1:
          status_line();
          break;
        case SIGUSR2:
          status_file();
          break;
      }

    if (istep < NTSTEP) {
      if (NPRSTEP * napkin.itstep > NPSTEP * napkin.itrstep) {
        disp_drift(driftfp);

        ++napkin.itrstep;
      }

      ++napkin.itstep;
    } else {
      est_gather_tde();

      // Bad.
      est_gather_x();

      if (NPRSTEP * napkin.ipstep > NPSTEP * napkin.iprstep) {
        disp_drift(driftfp);
        disp_tde(tdefp);

        ++napkin.iprstep;
      }

      ++napkin.ipstep;
    }
  }

  // Bad.
  printf("should be zero: %f +- %f\n", urgh.mean, urgh.std);

  if (fclose(tdefp) == EOF)
    err_abort(fclose);

  if (fclose(driftfp) == EOF)
    err_abort(fclose);

  status_line();
  status_file();

  gsl_rng_free(napkin.rng);
}

// TODO Split data files by dimension (use `sprintf`).
// TODO Abstract mean/var.
// TODO Fix V/K scaling and offsets.
// TODO Make periodicity conditional.
// TODO Abstract generic calculations for qho and He-4.
// TODO Figure out the correlation length for `stderr`.
// TODO Find home for lost souls.
// TODO Consider different masses for different polymers.
// TODO Try other Marsaglia's rngs (tried {U,V}NI with no observable effect).

int main(void) {
  not_main();

  return EXIT_SUCCESS;
}
