#include "floating.h"
#include "report.h"
#include "size.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define NDIM 2
#define NPOLY 32
#define NBEAD 42
#define NSTEP 12

static double const hbar = 1;

struct bead {
  double d[NDIM];
};

struct poly {
  struct bead r[NBEAD];
};

struct napkin {
  gsl_rng* rng;
  double m;
  double L;
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

*/

/*
Random initial configuration.
*/
static void randconf(void) {
  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead)
      for (size_t idim = 0;
          idim < NDIM;
          ++idim)
        napkin.R[ipoly].r[ibead].d[idim] =
          gsl_ran_flat(napkin.rng, 0, napkin.L);
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
      ++ipoly)
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

/*
Circular lattice initial configuration.
*/
static void clatconf(void) {
  double const w = ceil(pow(NPOLY, 1.0 / NDIM));
  double const v = napkin.L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
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

double sample(gsl_rng* const rng, double* const xs, size_t const n) {
  double x;
  gsl_ran_sample(rng, &x, 1, xs, n, sizeof x);
  return x;
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

static bool qho(void) {
  /* p. 44 */

  gsl_rng_env_setup();

  napkin.L = 5;

  napkin.m = 1;

  napkin.rng = gsl_rng_alloc(gsl_rng_mt19937);

  napkin.lambda = gsl_pow_2(hbar) / (2 * napkin.m);

  napkin.tau = NAN; // beta / M

  // randconf();
  // rlatconf();
  clatconf();

  FILE* const ensemfp = fopen("qho-ensemble.data", "w");
  if (ensemfp == NULL)
    halt(fopen);

  FILE* const napkinfp = fopen("qho-napkin.data", "w");
  if (napkinfp == NULL)
    halt(fopen);

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead) {
      fprintf(ensemfp, "%zu", ibead);

      for (size_t idim = 0;
          idim < NDIM;
          ++idim)
        fprintf(ensemfp, " %f", napkin.R[ipoly].r[ibead].d[idim]);

      fprintf(ensemfp, "\n");
    }

    fprintf(ensemfp, "\n");
  }

  fclose(napkinfp);

  fclose(ensemfp);

  gsl_rng_free(napkin.rng);

  return true;
}

int main(void) {
  reset();
  warn(reset);
  return qho() ? EXIT_SUCCESS : EXIT_FAILURE;
}
