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

#define NDIM ((size_t) 2)
#define NPOLY ((size_t) 4)
#define NBEAD ((size_t) 8)
#define NTSTEP ((size_t) 1 << 4)
#define NPSTEP ((size_t) 1 << 8)
#define NRSTEP ((size_t) 1 << 8)

/*
These assertions guarantee that adding or multiplying two steps or indices
never wraps around, which makes their manipulation easier.
*/
static_assert(NTSTEP <= SQRT_SIZE_MAX, "too many thermalization steps");
static_assert(NPSTEP <= SQRT_SIZE_MAX, "too many production steps");
static_assert(NRSTEP <= NPSTEP, "too many recording steps");

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
  size_t irstep; // Recording steps taken.
  double L; // Lattice constant.
  double beta;
  double lambda;
  double tau;
  double (* V2)(double); // Potential energy.
  double (* K2)(double); // Kinetic energy.
  double (* V2end)(double); // Potential energy at open ends.
  double (* K2end)(double); // Kinetic energy at open ends.
  struct {
    size_t N; // Samples.
    double mean; // Mean.
    double M2; // Second moment.
    double var; // Variance.
    double stderr; // Standard error of the mean.
  } tdE; // Thermodynamic estimator.
  struct poly R[NPOLY];
};

// Lost Souls () --------------------------------------------------------------

static size_t ran_index(gsl_rng* const rng, size_t const n) {
  return (size_t) gsl_rng_uniform_int(rng, n);
}

// Global State () ------------------------------------------------------------

static struct napkin napkin;

// Vector Math (bead) ---------------------------------------------------------

// This only takes the closest periodic image into account.
static struct bead bead_sub(struct bead const r0, struct bead const r1) {
  struct bead r;

  for (size_t idim = 0; idim < NDIM; ++idim)
    r.d[idim] = fp_periodic(r1.d[idim] - r0.d[idim], napkin.L);

  return r;
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

/*
Print parameters into stream `fp`.
*/
static void disp_drift(FILE* const fp) {
  fprintf(fp, "%zu %f %f\n",
      napkin.itstep + napkin.ipstep,
      napkin.params.ss.dx, napkin.params.comd.dx);
}

// Saving into Data Files (save)

static void save_esl(void) {
  FILE* const fp = fopen("qho-ensemble.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  disp_poly(fp);

  fclose(fp);
}

// Physical Moves (move) ------------------------------------------------------

static void move_accept_ss(void) {
  ++napkin.params.ss.accepted;
}

static void move_reject_ss(void) {
  napkin.R[napkin.history.ss.ipoly].r[napkin.history.ss.ibead] =
    napkin.history.ss.r;

  ++napkin.params.ss.rejected;
}

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
    napkin.R[ipoly].r[ibead].d[idim] = fp_wrap(
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
      napkin.R[ipoly].r[ibead].d[idim] = fp_wrap(
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
static double const epsilon = 1; // L--J 6--12 P
static double const sigma = 1; // L--J 6--12 P
static double const omega = 1; // HP

// System-Specific Quantities () ----------------------------------------------

static double Kbead(size_t const ipoly, size_t const ibead) {
  double K = 0;

  for (int i = 0; i < 2; ++i) {
    bool const p = ibead == 0 && i == 0;
    bool const q = ibead == NBEAD - 1 && i == 1;
    double (* const f)(double) = p || q ? napkin.K2end : napkin.K2;
    size_t const jbead = size_wrap(i == 0 ? ibead - 1 : ibead + 1, NBEAD);
    size_t const jpoly = p ? napkin.R[ipoly].from :
      q ? napkin.R[ipoly].to :
      ipoly;

    if (jpoly != SIZE_MAX)
      K += f(bead_norm2(bead_sub(napkin.R[ipoly].r[ibead],
              napkin.R[jpoly].r[jbead])));
  }

  return K;
}

static double Kpoly(size_t const ipoly) {
  double K = 0;

  for (size_t ibead = 0; ibead < NBEAD; ++ibead)
    K += Kbead(ipoly, ibead);

  return K;
}

static double Ktotal(void) {
  double K = 0;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    K += Kpoly(ipoly);

  return K;
}

static double Vbeads(size_t const ibead) {
  double V = 0;

  bool const p = ibead == 0 || ibead == NBEAD - 1;
  double (* const f)(double) = p ? napkin.V2end : napkin.V2;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < NPOLY; ++jpoly)
      V += f(bead_norm2(bead_sub(napkin.R[ipoly].r[ibead],
              napkin.R[jpoly].r[ibead])));

  return V;
}

static double Vtotal(void) {
  double V = 0;

  for (size_t ibead = 0; ibead < NBEAD; ++ibead)
    V += Vbeads(ibead);

  return V;
}

static double Stotal(void) {
  return Vtotal() + Ktotal();
}

// Estimators (est) -----------------------------------------------------------

static double est_tdE(void) {
  double const K = Ktotal() / (4 * napkin.lambda * napkin.tau);
  double const V = Vtotal() * napkin.tau;

  if (napkin.tdE.N == 0)
    return 0; // What is this heresy?

  return ((double) NDIM / 2 - K - V) / (napkin.tdE.N * napkin.tau);
}

static void est_gather_tdE(void) {
   double const E = est_tdE();
   ++napkin.tdE.N;
   double const Delta = E - napkin.tdE.mean;
   napkin.tdE.mean += Delta / (double) napkin.tdE.N;
   napkin.tdE.M2 += Delta * (E - napkin.tdE.mean);

   napkin.tdE.var = napkin.tdE.M2 / (double) (napkin.tdE.N - 1);
   napkin.tdE.stderr = sqrt(napkin.tdE.var / (double) napkin.tdE.N);
}

static double est_K(void) {
  double N = NAN, Rip1 = NAN, Ri = NAN; // ??
  double const cs = 2 * napkin.lambda * napkin.tau;
  return (NDIM * N * log(M_2PI * cs) + gsl_pow_2(Rip1 - Ri) / cs) / 2;
}

static double est_V(void) {
  double VRip1 = NAN, VRi = NAN; // ??
  return napkin.tau * (VRi + VRip1) / 2;
}

// Workers (work) -------------------------------------------------------------

static double work_ss(void) {
  size_t const ipoly = ran_index(napkin.rng, NPOLY);
  size_t const ibead = ran_index(napkin.rng, NBEAD);

  // Local `Vbeads`, because the rest stays the same.
  double const V0 = Vbeads(ibead);
  double const K0 = Kbead(ipoly, ibead);

  move_ss(ipoly, ibead);

  double const V1 = Vbeads(ibead);
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

// Other Crap () --------------------------------------------------------------

static void status_line(void) {
  printf("i / N = (%zu + %zu) / (%zu + %zu) = %zu %%\n"
      "E = %f +- %f\n",
      napkin.itstep, napkin.ipstep, NTSTEP, NPSTEP,
      100 * (napkin.itstep + napkin.ipstep) / (NTSTEP + NPSTEP),
      napkin.tdE.mean, napkin.tdE.stderr);
}

static void status_file(void) {
  save_esl();
}

static double lj612(double const r) {
  double const sigmar6 = gsl_pow_6(sigma / r);

  return 4 * epsilon * (gsl_pow_2(sigmar6) - sigmar6);
}

static double lj6122(double const r2) {
  double const sigmar2 = gsl_pow_2(sigma) / r2;

  return 4 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

static double harm(double const r) {
  return m * gsl_pow_2(omega * r) / 2;
}

static double harm2(double const r2) {
  return m * gsl_pow_2(omega) * r2 / 2;
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

  napkin.L = 5;
  napkin.beta = 1;
  napkin.lambda = gsl_pow_2(hbar) / (2 * m);
  napkin.tau = napkin.beta / NBEAD;

  napkin.V2 = lj6122;
  napkin.V2end = fp_zero;
  napkin.K2 = fp_identity;
  napkin.K2end = fp_identity;

  napkin.tdE.N = 0;
  napkin.tdE.mean = 0;
  napkin.tdE.M2 = 0;
  napkin.tdE.var = NAN; // This should go elsewhere.
  napkin.tdE.stderr = NAN; // This should go elsewhere.

  gsl_rng_env_setup();
  napkin.rng = gsl_rng_alloc(gsl_rng_default);

  conf_circlatt();
  force_close();

  FILE* const driftfp = fopen("qho-drift.data", "w");
  if (driftfp == NULL)
    err_abort(fopen);

  // Just to prevent empty data files...
  disp_drift(driftfp);
  fflush(driftfp);

  napkin.itstep = 0;
  napkin.ipstep = 0;
  napkin.irstep = 0;
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

    if (istep < NTSTEP)
      ++napkin.itstep;
    else {
      est_gather_tdE();

      if (NRSTEP * napkin.ipstep > NPSTEP * napkin.irstep) {
        disp_drift(driftfp);

        ++napkin.irstep;
      }

      ++napkin.ipstep;
    }
  }

  fclose(driftfp);

  status_line();
  status_file();

  gsl_rng_free(napkin.rng);
}

// abstract mean/var, fix V/K scaling, find home for lost souls

int main(void) {
  not_main();

  return EXIT_SUCCESS;
}
