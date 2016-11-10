#include "floating.h"
#include "report.h"
#include "size.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define NDIM 2
#define NPOLY 3
#define NBEAD 4
#define NSTEP 12

struct pos {
  double coords[NDIM];
};

struct poly {
  struct pos beads[NBEAD];
};

struct ensemble {
  gsl_rng* rng;
  double lambda;
  double tau;
  struct poly polys[NPOLY];
};

struct napkin {
  double energy;
};

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

static void bisectmv(struct ensemble* const ens,
    size_t const ipoly, size_t const ibead) {
  for (size_t idim = 0;
      idim < NDIM;
      ++idim) {
    /*
    r = gsl_ran_gaussian(rng, sqrt(ens.lambda * ens.tau)) +
      midpoint(ens.polys[ipoly].beads[size_wrap(ibead - 1)].coords[idim],
          ens.polys[ipoly].beads[size_wrap(ibead + 1)].coords[idim]);
    ens.polys[ipoly].beads[ibead].coords[idim] = ...;
    */
  }
}

static void mlevmv(struct ensemble* const ens,
    size_t const ipoly, size_t const ibead, size_t const nlev) {
}

static void displmv(struct ensemble* const ens, size_t const ipoly) {
}

static double Kthing(struct ensemble* const ens,
    size_t const ipoly, size_t const ibead) {
  double const cse = 2 * lambda * tau;

  return -(d * N / 2) * log(M_2PI * cse)
    + gsl_pow_2(length(R[m - 1] - R[m])) / (2 * cse);
}

// prim approx
static double Uthing(struct ensemble* const ens,
    size_t const ipoly, size_t const ibead) {
  return (tau / 2) * (V(R[m - 1]) + V(R[m]));
}

// prim approx also
static double Sthing(struct ensemble* const ens,
    size_t const ipoly, size_t const ibead) {
  return Kthing(ens, ipoly, ibead) + Uthing(ens, ipoly, ibead);
}

bool qho(void) {
  struct ensemble ens;
  struct napkin nk;

  /* p. 44 */

  gsl_rng_env_setup();

  ens.rng = gsl_rng_alloc(gsl_rng_mt19937);

  ens.lambda = NAN; // hbar^2 / 2 m

  ens.tau = NAN; // beta / M

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly)
    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead)
      for (size_t idim = 0;
          idim < NDIM;
          ++idim)
        ens.polys[ipoly].beads[ibead].coords[idim] = gsl_rng_uniform(rng);

  FILE* const ensfp = fopen("qho-ensemble.data", "w");
  if (ensfp == NULL)
    halt(fopen);

  FILE* const napfp = fopen("qho-napkin.data", "w");
  if (napfp == NULL)
    halt(fopen);

  for (size_t ipoly = 0;
      ipoly < NPOLY;
      ++ipoly) {
    for (size_t ibead = 0;
        ibead < NBEAD;
        ++ibead) {
      struct pos xy = ens.polys[ipoly].beads[ibead];

      fprintf(ensfp, "%zu %f %f\n", ibead, xy.coords[0], xy.coords[1]);
    }

    fprintf(ensfp, "\n");
  }

  fclose(napfp);

  fclose(ensfp);

  gsl_rng_free(rng);

  return true;
}
