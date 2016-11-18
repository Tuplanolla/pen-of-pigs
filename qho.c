#include "floating.h"
#include "report.h"
#include "size.h"
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
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
#define NPSTEP 256
#define NRSTEP 100

static_assert(NRSTEP <= NPSTEP,
    "too many recording steps");

// This is just a static version of `assert(NPSTEP <= sqrt(SIZE_MAX))`.
static_assert(NPSTEP <= (size_t) 1 << CHAR_BIT * sizeof (size_t) / 2,
    "too many production steps");

// Constants () ---------------------------------------------------------------

static double const hbar = 1;

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
  struct poly R[NPOLY];
};

struct paper {
  double energy;
};

// Global State () ------------------------------------------------------------

static struct napkin napkin;

static struct paper paper;

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
static void dispdisp(FILE* const fp) {
  fprintf(fp, "%f\n", napkin.dx);
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

// Physical Adjustments (adjust) ----------------------------------------------

static void adjust_dx(void) {
  napkin.dx *= (double) (napkin.accepted + 1) / (double) (napkin.rejected + 1);
}

// System-Specific Quantities () ----------------------------------------------

static double Kone(size_t const ipoly, size_t const ibead) {
  /*
  double const cse = 2 * lambda * tau;

  return -(d * N / 2) * log(M_2PI * cse)
    + gsl_pow_2(length(R[m - 1] - R[m])) / (2 * cse);
  */
  return 0;
}

static double Vone(size_t const ipoly, size_t const ibead) {
  /*
  return (tau / 2) * (V(R[m - 1]) + V(R[m]));
  */
  return 0;
}

static double Sone(size_t const ipoly, size_t const ibead) {
  /*
  return Kone(napkin, ipoly, ibead) + Uone(napkin, ipoly, ibead);
  */
  return 0;
}

static double Kall(void) {
  return 0;
}

static double Vall(void) {
  return 0;
}

static double Sall(void) {
  return 0;
}

// Work Horses () -------------------------------------------------------------

static void subwork(void) {
  move_trial((size_t) gsl_rng_uniform_int(napkin.rng, NPOLY),
      (size_t) gsl_rng_uniform_int(napkin.rng, NBEAD));
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

static void work(void) {
  double const T0 = 0.5;
  double const K0 = 0.5;

  subwork();

  double const T1 = 0.5;
  double const K1 = 0.5;

  double const DeltaS_T = napkin.tau * (T1 - T0);
  double const DeltaS_K = (1 / (4 * napkin.lambda * napkin.tau)) * (K1 - K0);

  choice(DeltaS_T + DeltaS_K);

  // adjust_dx();
}

static void nop(void) {}

int main(void) {
  /* p. 44 */

  reset();

  gsl_rng_env_setup();

  napkin.L = 5;
  napkin.m = 1;
  napkin.dx = 1;
  napkin.beta = 1;

  napkin.rng = gsl_rng_alloc(gsl_rng_mt19937);

  napkin.lambda = gsl_pow_2(hbar) / (2 * napkin.m);

  napkin.tau = napkin.beta / NBEAD;

  napkin.undo = nop;

  conf_circlatt();

  FILE* const dispfp = fopen("qho-displacement.data", "w");
  if (dispfp == NULL)
    halt(fopen);

  // Thermalize.
  for (size_t itstep = 0;
      itstep < NTSTEP;
      ++itstep)
    work();

  // Work productively.
  for (size_t ipstep = 0, irstep = 0;
      ipstep < NPSTEP;
      ++ipstep) {
    work();

    // Update paper with estimators

    if (NRSTEP * ipstep > NPSTEP * irstep) {
      // Save results into files
      dispdisp(dispfp);

      ++irstep;
    }
  }

  // Closing index exchange test...
  napkin.R[0].i = 2;
  napkin.R[2].i = 0;

  fclose(dispfp);

  FILE* const fp = fopen("qho-ensemble.data", "w");
  if (fp == NULL)
    halt(fopen);

  disp_poly(fp);

  fclose(fp);

  gsl_rng_free(napkin.rng);

  return EXIT_SUCCESS;
}
