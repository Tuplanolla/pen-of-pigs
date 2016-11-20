#include "apxtime.h"
#include "err.h"
#include "flt.h"
#include "phys.h"
#include "sigs.h"
#include "size.h"
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <signal.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Static Constants (N) -------------------------------------------------------

#define NDIM 2
#define NPOLY 4
#define NBEAD 8
#define NTSTEP (1 << 4)
#define NPSTEP (1 << 16)
#define NRSTEP (1 << 8)

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
    double d[NDIM];
  } trial;
  struct {
    int dummy;
  } bisect;
};

struct napkin {
  gsl_rng* rng;
  size_t accepted;
  size_t rejected;
  void (* undo)(void);
  union hist history;
  double dx;
  double m;
  double L;
  double beta;
  double lambda;
  double tau;
  double (* V2)(double);
  double (* V2cl)(double);
  double (* K2)(double);
  double (* K2cl)(double);
  double energy;
  struct poly R[NPOLY];
};

// Global State () ------------------------------------------------------------

static struct napkin napkin;

// Vector Math (bead) ---------------------------------------------------------

static struct bead bead_add(struct bead const r0, struct bead const r1) {
  struct bead r;

  for (size_t idim = 0; idim < NDIM; ++idim)
    r.d[idim] = r0.d[idim] + r1.d[idim];

  return r;
}

static struct bead bead_sub(struct bead const r0, struct bead const r1) {
  struct bead r;

  for (size_t idim = 0; idim < NDIM; ++idim)
    r.d[idim] = r0.d[idim] - r1.d[idim];

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
  size_t const ipoly = napkin.history.trial.ipoly;
  size_t const ibead = napkin.history.trial.ibead;

  for (size_t idim = 0; idim < NDIM; ++idim)
    napkin.R[ipoly].r[ibead].d[idim] = napkin.history.trial.d[idim];
}

static void move_trial(size_t const ipoly, size_t const ibead) {
  napkin.undo = unmove_trial;

  napkin.history.trial.ipoly = ipoly;
  napkin.history.trial.ibead = ibead;

  for (size_t idim = 0; idim < NDIM; ++idim) {
    napkin.history.trial.d[idim] = napkin.R[ipoly].r[ibead].d[idim];

    napkin.R[ipoly].r[ibead].d[idim] = wrap(
        napkin.R[ipoly].r[ibead].d[idim] +
        napkin.dx * gsl_ran_flat(napkin.rng, -1, 1),
      napkin.L);
  }
}

static void move_bisect(size_t const ipoly, size_t const ibead) {
  for (size_t idim = 0; idim < NDIM; ++idim) {
    napkin.R[ipoly].r[ibead].d[idim] =
      gsl_ran_gaussian(napkin.rng, sqrt(napkin.lambda * napkin.tau)) +
      midpoint(napkin.R[ipoly].r[size_wrap(ibead - 1, NBEAD)].d[idim],
          napkin.R[ipoly].r[size_wrap(ibead + 1, NBEAD)].d[idim]);
  }
}

static void move_multlev(size_t const ipoly, size_t const ibead,
    size_t const nlev) {
}

static void move_displace(size_t const ipoly) {
  for (size_t idim = 0; idim < NDIM; ++idim) {
    double const r =
      gsl_ran_gaussian(napkin.rng, sqrt(napkin.lambda * napkin.tau));

    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      napkin.R[ipoly].r[ibead].d[idim] += r;
  }
}

// Adjustments (adjust) -------------------------------------------------------

static double ratio(void) {
  return (double) (napkin.accepted + 1) / (double) (napkin.rejected + 1);
}

static void adjust_always(void) {
  napkin.dx *= ratio();
}

static void adjust_factor(double const factor) {
  double const p = ratio();

  napkin.dx *= p > factor || p < 1 / factor ? p : 1;
}

// Constants () ---------------------------------------------------------------

static double const hbar = 1;
static double const epsilon = 1;
static double const sigma = 1;

// System-Specific Quantities () ----------------------------------------------

static double Sone(size_t const ipoly, size_t const ibead) {
  /*
  return Kone(napkin, ipoly, ibead) + Uone(napkin, ipoly, ibead);
  */
  return 0;
}

static double Stotal(void) {
  return 0;
}

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

static double Ktotal(void) {
  double K = 0;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      K += Kbead(ipoly, ibead);

  return K;
}

static double Vbead(size_t const ibead) {
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
    V += Vbead(ibead);

  return V;
}

// Work Horses () -------------------------------------------------------------

static double subwork(void) {
  size_t const ipoly = (size_t) gsl_rng_uniform_int(napkin.rng, NPOLY);
  size_t const ibead = (size_t) gsl_rng_uniform_int(napkin.rng, NBEAD);

  // Local `Vbead`, because the rest stays the same.
  double const T0 = Vbead(ibead);
  double const K0 = Kbead(ipoly, ibead);

  move_trial(ipoly, ibead);

  double const T1 = Vbead(ibead);
  double const K1 = Kbead(ipoly, ibead);

  double const DeltaS_T = napkin.tau * (T1 - T0);
  double const DeltaS_K = (1 / (4 * napkin.lambda * napkin.tau)) * (K1 - K0);

  return DeltaS_T + DeltaS_K;
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

static void status(void) {
  save_esl();
}

static void work(void) {
  choice(subwork());

  adjust_factor(1.25);

  int signum;
  if (sigs_use(&signum) && signum == SIGUSR1)
    status();
}

static double Vee2(double const r2) {
  return lj6122(r2, epsilon, gsl_pow_2(sigma));
}

static void not_main(void) {
  err_reset();

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

  napkin.L = 5;
  napkin.m = 1;
  napkin.dx = 2;
  napkin.beta = 1;
  napkin.V2 = Vee2;
  napkin.V2cl = zero;
  napkin.K2 = identity;
  napkin.K2cl = identity;

  gsl_rng_env_setup();
  napkin.rng = gsl_rng_alloc(gsl_rng_default);

  napkin.lambda = gsl_pow_2(hbar) / (2 * napkin.m);

  napkin.tau = napkin.beta / NBEAD;

  napkin.undo = unmove_null;

  conf_circlatt();
  force_close();

  FILE* const driftfp = fopen("qho-drift.data", "w");
  if (driftfp == NULL)
    err_abort(fopen);

  disp_drift(driftfp); // Just to prevent empty data files.

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

  status();

  gsl_rng_free(napkin.rng);
}

int main(void) {
  not_main();

  return EXIT_SUCCESS;
}
