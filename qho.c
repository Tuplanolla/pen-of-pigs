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
#define NSTAB ((size_t) 1 << 6)
#define NTSTEP ((size_t) 1 << 4)
#define NPSTEP ((size_t) 1 << 16)
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

union hist {
  struct {
    size_t ipoly;
    size_t ibead;
    struct bead r;
  } trial;
  struct {
    size_t ipoly;
    struct poly R;
  } displace;
};

struct napkin {
  gsl_rng* rng;
  size_t accepted; // These should be per-move.
  size_t rejected; // Yes.
  double dx; // Indeed.
  void (* undo)(void);
  union hist history;
  double L;
  double beta;
  double lambda;
  double tau;
  double (* V2)(double);
  double (* K2)(double);
  double (* V2cl)(double); // Closure properties.
  double (* K2cl)(double); // These may be bad for assumptions in displacement.
  double energy;
  struct poly R[NPOLY];
};

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
Print trial displacement into stream `fp`.
*/
static void disp_drift(FILE* const fp) {
  fprintf(fp, "%zu %f\n", napkin.accepted + napkin.rejected, napkin.dx);
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

static void unmove_null(void) {}

static void move_null(void) {}

static void unmove_trial(void) {
  napkin.R[napkin.history.trial.ipoly].r[napkin.history.trial.ibead] =
    napkin.history.trial.r;
}

static void move_trial(size_t const ipoly, size_t const ibead) {
  napkin.undo = unmove_trial;

  napkin.history.trial.ipoly = ipoly;
  napkin.history.trial.ibead = ibead;
  napkin.history.trial.r = napkin.R[ipoly].r[ibead];

  for (size_t idim = 0; idim < NDIM; ++idim)
    napkin.R[ipoly].r[ibead].d[idim] = fp_wrap(
        napkin.R[ipoly].r[ibead].d[idim] +
        napkin.dx * gsl_ran_flat(napkin.rng, -1, 1),
        napkin.L);
}

static void unmove_displace(void) {
  napkin.R[napkin.history.displace.ipoly] = napkin.history.displace.R;
}

static void move_displace(size_t const ipoly) {
  napkin.undo = unmove_displace;

  napkin.history.displace.ipoly = ipoly;
  napkin.history.displace.R = napkin.R[ipoly];

  for (size_t idim = 0; idim < NDIM; ++idim) {
    double const x = napkin.dx * gsl_ran_flat(napkin.rng, -1, 1);

    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      napkin.R[ipoly].r[ibead].d[idim] = fp_wrap(
          napkin.R[ipoly].r[ibead].d[idim] + x,
          napkin.L);
  }
}

static void unmove_bisect(size_t const ipoly, size_t const ibead) {
  err_abort(NULL);
}

static void move_bisect(size_t const ipoly, size_t const ibead) {
  err_abort(NULL);

  for (size_t idim = 0; idim < NDIM; ++idim) {
    napkin.R[ipoly].r[ibead].d[idim] =
      gsl_ran_gaussian(napkin.rng, sqrt(napkin.lambda * napkin.tau)) +
      fp_midpoint(napkin.R[ipoly].r[size_wrap(ibead - 1, NBEAD)].d[idim],
          napkin.R[ipoly].r[size_wrap(ibead + 1, NBEAD)].d[idim]);
  }
}

// Adjustments (adjust) -------------------------------------------------------

static double ratio(void) {
  return (double) (napkin.accepted + NSTAB) /
    (double) (napkin.rejected + NSTAB);
}

static size_t total(void) {
  return napkin.accepted + napkin.rejected;
}

static void adjust_yes(void) {
  napkin.dx = fp_wrap(napkin.dx * ratio(), napkin.L);
}

static void adjust_maybe(void) {
  if (total() % NSTAB == 0)
    adjust_yes();
}

// Constants () ---------------------------------------------------------------

static double const hbar = 1;
static double const m = 1;
static double const epsilon = 1; // L--J 6--12 P
static double const sigma = 1; // L--J 6--12 P
static double const omega = 1; // HP

// Estimators (est) -----------------------------------------------------------

static double est_K(void) {
  double N = NAN, Rip1 = NAN, Ri = NAN; // ??
  double const cs = 2 * napkin.lambda * napkin.tau;
  return (NDIM * N * log(M_2PI * cs) + gsl_pow_2(Rip1 - Ri) / cs) / 2;
}

static double est_V(void) {
  double VRip1 = NAN, VRi = NAN; // ??
  return napkin.tau * (VRi + VRip1) / 2;
}

// System-Specific Quantities () ----------------------------------------------

static double Kbead(size_t const ipoly, size_t const ibead) {
  double K = 0;

  for (int i = 0; i < 2; ++i) {
    bool const p = ibead == 0 && i == 0;
    bool const q = ibead == NBEAD - 1 && i == 1;
    double (* const f)(double) = p || q ? napkin.K2cl : napkin.K2;
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
  double (* const f)(double) = p ? napkin.V2cl : napkin.V2;

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

// Work Horses () -------------------------------------------------------------

typedef double (* worker)(void);

static double work_trial(void) {
  size_t const ipoly = (size_t) gsl_rng_uniform_int(napkin.rng, NPOLY);
  size_t const ibead = (size_t) gsl_rng_uniform_int(napkin.rng, NBEAD);

  // Local `Vbeads`, because the rest stays the same.
  double const T0 = Vbeads(ibead);
  double const K0 = Kbead(ipoly, ibead);

  move_trial(ipoly, ibead);

  double const T1 = Vbeads(ibead);
  double const K1 = Kbead(ipoly, ibead);

  double const DeltaS_T = napkin.tau * (T1 - T0);
  double const DeltaS_K = (1 / (4 * napkin.lambda * napkin.tau)) * (K1 - K0);

  return DeltaS_T + DeltaS_K;
}

static double work_displace(void) {
  size_t const ipoly = (size_t) gsl_rng_uniform_int(napkin.rng, NPOLY);

  // Global `Vtotal`, because `K` is unchanged (maybe not for permutations).
  double const T0 = Vtotal();

  move_displace(ipoly);

  double const T1 = Vtotal();

  double const DeltaS_T = napkin.tau * (T1 - T0);

  return DeltaS_T;
}

static void choice(double const DeltaS) {
  if (DeltaS <= 0)
    ++napkin.accepted;
  else {
    double const q = gsl_ran_flat(napkin.rng, 0, 1);

    if (q <= exp(-DeltaS))
      ++napkin.accepted;
    else {
      napkin.undo();

      ++napkin.rejected;
    }
  }
}

static void status_line(void) {
  printf("z = (%zu + %zu) / (%zu + %zu) = %zu %%\ndx = %f\n",
      napkin.accepted, napkin.rejected, NTSTEP, NPSTEP,
      100 * (napkin.accepted + napkin.rejected) / (NTSTEP + NPSTEP),
      napkin.dx);
}

static void status_file(void) {
  save_esl();
}

static void work(void) {
  worker const workers[] = {work_trial, work_displace};

  size_t const iworker = (size_t) gsl_rng_uniform_int(napkin.rng,
      sizeof workers / sizeof *workers);
  choice(workers[iworker]());

  adjust_maybe();

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

static void not_main(void) {
  err_reset();

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

  napkin.accepted = 0;
  napkin.rejected = 0;
  napkin.dx = 1;
  napkin.L = 5;
  napkin.beta = 1;
  napkin.lambda = gsl_pow_2(hbar) / (2 * m);
  napkin.tau = napkin.beta / NBEAD;
  napkin.V2 = lj6122;
  napkin.V2cl = fp_zero;
  napkin.K2 = fp_identity;
  napkin.K2cl = fp_identity;

  gsl_rng_env_setup();
  napkin.rng = gsl_rng_alloc(gsl_rng_default);

  napkin.undo = unmove_null;

  conf_circlatt();
  force_close();

  FILE* const driftfp = fopen("qho-drift.data", "w");
  if (driftfp == NULL)
    err_abort(fopen);

  // Just to prevent empty data files...
  disp_drift(driftfp);
  fflush(driftfp);

  for (size_t istep = 0, itstep = 0, ipstep = 0, irstep = 0;
      istep < NTSTEP + NPSTEP;
      ++istep) {
    work();

    if (istep < NTSTEP)
      ++itstep;
    else {
      // Update napkin with estimators

      if (NRSTEP * ipstep > NPSTEP * irstep) {
        // Save results into files
        disp_drift(driftfp);

        ++irstep;
      }

      ++ipstep;
    }
  }

  fclose(driftfp);

  status_file();

  gsl_rng_free(napkin.rng);
}

int main(void) {
  not_main();

  return EXIT_SUCCESS;
}
