#include "apxtime.h"
#include "err.h"
#include "flt.h"
#include "sig.h"
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

#define NDIM 3
#define NPOLY 4
#define NBEAD 8
#define NTSTEP 16
#define NPSTEP (1 << 20)
#define NRSTEP 100

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
  size_t i;
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
  double dt;
  double t0;
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
  double energy;
  struct poly R[NPOLY];
};

// Global State () ------------------------------------------------------------

static struct napkin napkin;

// Vector Math (??) -----------------------------------------------------------

static void add(double* const d,
    double const* const d0, double const* const d1) {
  for (size_t idim = 0;
      idim < NDIM;
      ++idim)
    d[idim] = d0[idim] + d1[idim];
}

static void sub(double* const d,
    double const* const d0, double const* const d1) {
  for (size_t idim = 0;
      idim < NDIM;
      ++idim)
    d[idim] = d1[idim] - d0[idim];
}

static double normsq(double const* const d) {
  double s = 0;

  for (size_t idim = 0;
      idim < NDIM;
      ++idim)
    s += gsl_pow_2(d[idim]);

  return s;
}

static double norm(double const* const d) {
  return sqrt(normsq(d));
}

// Initial Configurations (conf) ----------------------------------------------

/*
Random initial configuration.
*/
static void conf_rand(void) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    napkin.R[ipoly].i = ipoly;

    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead)
      for (size_t idim = 0;
          idim < NDIM;
          ++idim)
        napkin.R[ipoly].r[ibead].d[idim] =
          gsl_ran_flat(napkin.rng, 0, napkin.L);
  }
}

/*
Random point initial configuration.
*/
static void conf_randpt(void) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    napkin.R[ipoly].i = ipoly;

    for (size_t idim = 0;
        idim < NDIM;
        ++idim) {
      double const x = gsl_ran_flat(napkin.rng, 0, napkin.L);

      for (size_t ibead = 0;
          ibead < NBEAD;
          ++ibead)
        napkin.R[ipoly].r[ibead].d[idim] = x;
    }
  }
}

/*
Random lattice initial configuration.
*/
static void conf_randlatt(void) {
  double const w = ceil(pow(NPOLY, 1.0 / NDIM));
  double const v = napkin.L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    napkin.R[ipoly].i = ipoly;

    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead) {
      size_t m = ipoly;

      for (size_t idim = 0;
          idim < NDIM;
          ++idim) {
        size_div_t const z = size_div(m, n);

        napkin.R[ipoly].r[ibead].d[idim] =
          (double) z.rem * v + gsl_ran_flat(napkin.rng, 0, v);
        m = z.quot;
      }
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

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    napkin.R[ipoly].i = ipoly;

    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead) {
      size_t m = ipoly;

      for (size_t idim = 0;
          idim < NDIM;
          ++idim) {
        size_div_t const z = size_div(m, n);

        napkin.R[ipoly].r[ibead].d[idim] = (double) z.rem * v;
        m = z.quot;
      }
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

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    napkin.R[ipoly].i = ipoly;

    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead) {
      double const t = (double) ibead / NBEAD;

      size_t m = ipoly;

      for (size_t idim = 0;
          idim < NDIM;
          ++idim) {
        double const phi = (double) idim / NDIM;
        double const y = sin(M_2PI * (t + phi / 2)) / 2;

        size_div_t const z = size_div(m, n);
        napkin.R[ipoly].r[ibead].d[idim] = ((double) z.rem + (y + 1) / 2) * v;
        m = z.quot;
      }
    }
  }
}

// Configuration Forcing (force) ----------------------------------------------

/*
Close the current configuration.
*/
static void force_close(void) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
    napkin.R[ipoly].i = ipoly;
}

/*
Open the current configuration.
*/
static void force_open(void) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
    napkin.R[ipoly].i = SIZE_MAX;
}

// Printing into Data Files (disp) --------------------------------------------

/*
Print bead into stream `fp`.
*/
static void disp_bead(FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  fprintf(fp, "%zu", iindex);

  for (size_t idim = 0;
      idim < NDIM;
      ++idim)
    fprintf(fp, " %f", napkin.R[ipoly].r[ibead].d[idim]);

  fprintf(fp, "\n");
}

/*
Print polymer into stream `fp`.
*/
static void disp_poly(FILE* const fp) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead)
      disp_bead(fp, ibead, ipoly, ibead);

    if (napkin.R[ipoly].i != SIZE_MAX)
      disp_bead(fp, NBEAD, napkin.R[ipoly].i, 0);

    fprintf(fp, "\n");
  }
}

/*
Print trial displacement into stream `fp`.
*/
static void disp_drift(FILE* const fp) {
  fprintf(fp, "%f\n", napkin.dx);
}

// Saving into Data Files (save)

static void save_esl(void) {
  FILE* const fp = fopen("qho-ensemble.data", "w");
  if (fp == NULL)
    halt(fopen);

  disp_poly(fp);

  fclose(fp);
}

// Physical Moves (move) ------------------------------------------------------

static void unmove_trial(void) {
  size_t const ipoly = napkin.history.trial.ipoly;
  size_t const ibead = napkin.history.trial.ibead;

  for (size_t idim = 0;
      idim < NDIM;
      ++idim)
    napkin.R[ipoly].r[ibead].d[idim] = napkin.history.trial.d[idim];
}

static void move_trial(size_t const ipoly, size_t const ibead) {
  napkin.undo = unmove_trial;

  napkin.history.trial.ipoly = ipoly;
  napkin.history.trial.ibead = ibead;

  for (size_t idim = 0;
      idim < NDIM;
      ++idim) {
    napkin.history.trial.d[idim] = napkin.R[ipoly].r[ibead].d[idim];

    napkin.R[ipoly].r[ibead].d[idim] = wrap(
        napkin.R[ipoly].r[ibead].d[idim] +
        napkin.dx * gsl_ran_flat(napkin.rng, -1, 1),
      napkin.L);
  }
}

static void move_bisect(size_t const ipoly, size_t const ibead) {
  for (size_t idim = 0;
      idim < NDIM;
      ++idim) {
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
  for (size_t idim = 0;
      idim < NDIM;
      ++idim) {
    double const r =
      gsl_ran_gaussian(napkin.rng, sqrt(napkin.lambda * napkin.tau));

    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead)
      napkin.R[ipoly].r[ibead].d[idim] += r;
  }
}

// Adjustments (adjust) -------------------------------------------------------

static void adjust_dx(void) {
  napkin.dx *= (double) (napkin.accepted + 1) / (double) (napkin.rejected + 1);
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

static double lj612(double const r) {
  return 4 * epsilon *
    (gsl_pow_uint(sigma / r, 12) - gsl_pow_6(sigma / r));
}

static double lj6122(double const r2) {
  double const sigma2 = gsl_pow_2(sigma);

  return 4 * epsilon *
    (gsl_pow_6(sigma2 / r2) - gsl_pow_3(sigma2 / r2));
}

static double Kbead(size_t const ipoly, size_t const ibead) {
  double K = 0;

  // Note: self-closure assumed!
  size_t const xs[] = {1, NBEAD - 1};
  for (size_t i = 0;
      i < sizeof xs / sizeof *xs;
      ++i) {
    double d[NDIM];
    sub(d, napkin.R[ipoly].r[ibead].d,
        napkin.R[ipoly].r[size_wrap(ibead + xs[i], NBEAD)].d);

    K += normsq(d);
  }

  return K;
}

static double Ktotal(void) {
  double K = 0;

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead)
      K += Kbead(ipoly, ibead);

  return K;
}

static double Vbead(size_t const ibead) {
  double V = 0;

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
    for (size_t jpoly = ipoly + 1;
        jpoly < NPOLY;
        ++jpoly) {
      double d[NDIM];
      sub(d, napkin.R[ipoly].r[ibead].d, napkin.R[jpoly].r[ibead].d);

      V += lj612(norm(d));
    }

  return V;
}

static double Vtotal(void) {
  double V = 0;

  for (size_t ibead = 0;
      ibead < NBEAD;
      ++ibead)
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

  if ((napkin.accepted + napkin.rejected) % 256 == 0)
    adjust_dx();

  int signum;
  if (sig_use(&signum) && signum == SIGUSR1)
    status();
}

static void nop(void) {}

static void not_main(void) {
  reset();

  int const xs[] = {SIGUSR1, SIGUSR2};
  if (sig_register(xs, sizeof xs / sizeof *xs) != SIZE_MAX)
    halt(sig_register);

  napkin.dt = 1;
  napkin.t0 = now();

  napkin.L = 5;
  napkin.m = 1;
  napkin.dx = 2;
  napkin.beta = 1;

  gsl_rng_env_setup();
  napkin.rng = gsl_rng_alloc(gsl_rng_default);

  napkin.lambda = gsl_pow_2(hbar) / (2 * napkin.m);

  napkin.tau = napkin.beta / NBEAD;

  napkin.undo = nop;

  conf_circlatt();

  FILE* const driftfp = fopen("qho-drift.data", "w");
  if (driftfp == NULL)
    halt(fopen);

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

  save_esl();

  gsl_rng_free(napkin.rng);
}

int main(void) {
  not_main();

  return EXIT_SUCCESS;
}
