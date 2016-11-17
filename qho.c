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

#define NDIM 2
#define NPOLY 4
#define NBEAD 8
#define NTSTEP 16
#define NPSTEP 256
#define NRSTEP 100

static double const hbar = 1;

struct bead {
  double d[NDIM];
};

struct poly {
  size_t i;
  struct bead r[NBEAD];
};

struct napkin {
  gsl_rng* rng;
  size_t accepted;
  size_t rejected;
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

static struct napkin napkin;

static struct paper paper;

/*

double linkact() {
}

double action() {
  return -log(linkact());
}

double denmat() {
  return integral(exp(-sum(action())));
}

double sample(gsl_rng* const rng, double* const xs, size_t const n) {
  double x;
  gsl_ran_sample(rng, &x, 1, xs, n, sizeof x);
  return x;
}

*/

/*
Random initial configuration.
*/
static void rconf(void) {
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
static void rpconf(void) {
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
static void rlatconf(void) {
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
static void platconf(void) {
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
static void clatconf(void) {
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

/*
Close the current configuration.
*/
static void closec(void) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
    napkin.R[ipoly].i = ipoly;
}

/*
Open the current configuration.
*/
static void openc(void) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
    napkin.R[ipoly].i = SIZE_MAX;
}

/*
Print bead into stream `fp`.
*/
static void dispbead(FILE* const fp,
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
static void disppoly(FILE* const fp) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead)
      dispbead(fp, ibead, ipoly, ibead);

    if (napkin.R[ipoly].i != SIZE_MAX)
      dispbead(fp, NBEAD, napkin.R[ipoly].i, 0);

    fprintf(fp, "\n");
  }
}

/*
Print trial displacement into stream `fp`.
*/
static void dispdisp(FILE* const fp) {
  fprintf(fp, "%f\n", napkin.dx);
}

static void trialmv(size_t const ipoly, size_t const ibead) {
  for (size_t idim = 0;
      idim < NDIM;
      ++idim) {
    napkin.R[ipoly].r[ibead].d[idim] = wrap(
        napkin.R[ipoly].r[ibead].d[idim] +
        napkin.dx * gsl_ran_flat(napkin.rng, -1, 1),
      napkin.L);
  }
  // Displacement adjustment test...
  ++napkin.accepted;
  ++napkin.rejected;
  napkin.accepted += (size_t) gsl_rng_uniform_int(napkin.rng, 257);
  napkin.rejected += (size_t) gsl_rng_uniform_int(napkin.rng, 242);
}

static void adjustdx(void) {
  napkin.dx *= (double) (napkin.accepted + 1) / (double) (napkin.rejected + 1);
}

static void bisectmv(size_t const ipoly, size_t const ibead) {
  for (size_t idim = 0;
      idim < NDIM;
      ++idim) {
    napkin.R[ipoly].r[ibead].d[idim] =
      gsl_ran_gaussian(napkin.rng, sqrt(napkin.lambda * napkin.tau)) +
      midpoint(napkin.R[ipoly].r[size_wrap(ibead - 1, NBEAD)].d[idim],
          napkin.R[ipoly].r[size_wrap(ibead + 1, NBEAD)].d[idim]);
  }
}

static void mullevmv(size_t const ipoly, size_t const ibead, size_t const nlev) {
}

static void displmv(size_t const ipoly) {
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

/*

static double Kthing(size_t const ipoly, size_t const ibead) {
  double const cse = 2 * lambda * tau;

  return -(d * N / 2) * log(M_2PI * cse)
    + gsl_pow_2(length(R[m - 1] - R[m])) / (2 * cse);
}

// prim approx
static double Uthing(size_t const ipoly, size_t const ibead) {
  return (tau / 2) * (V(R[m - 1]) + V(R[m]));
}

// prim approx also
// but wrong
static double Sthing(size_t const ipoly, size_t const ibead) {
  return Kthing(napkin, ipoly, ibead) + Uthing(napkin, ipoly, ibead);
}

*/

static void work(void) {
  trialmv((size_t) gsl_rng_uniform_int(napkin.rng, NPOLY),
      (size_t) gsl_rng_uniform_int(napkin.rng, NBEAD));

  // Compute the diffence in potential action between the new and old config

  // Accept or reject

  adjustdx();
}

int main(void) {
  /* p. 44 */

  reset();

  assert(NRSTEP <= NPSTEP);
  assert(NPSTEP <= (size_t) 1 << CHAR_BIT * sizeof (size_t) / 2);

  gsl_rng_env_setup();

  napkin.L = 5;
  napkin.m = 1;
  napkin.dx = 1;
  napkin.beta = 1;

  napkin.rng = gsl_rng_alloc(gsl_rng_mt19937);

  napkin.lambda = gsl_pow_2(hbar) / (2 * napkin.m);

  napkin.tau = napkin.beta / NBEAD;

  clatconf();

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

  disppoly(fp);

  fclose(fp);

  gsl_rng_free(napkin.rng);

  return EXIT_SUCCESS;
}
