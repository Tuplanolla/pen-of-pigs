#include "floating.h"
#include "nth.h"
#include "report.h"
#include "size.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define NDIM 3
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

static double ljfreq(size_t const n) {
  return M_2PI * (double) (n < 2 ? 1 : nth_prime(n - 2));
}

static double ljphase(size_t const n) {
  return (M_2PI / 2) * ((double) n / NDIM);
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
        /*
        d[i] = sincos(n[i] * M_2PI * t + k[i] * M_PI)
        n[i] -- prime (rel-coprime)
        n[i] * k[j] - n[j] * k[i] -- not integer
        ----
        Example for 3.
        p[0] * k[1] - p[1] * k[0] -- not integer
        p[0] * k[2] - p[2] * k[0] -- not integer
        p[1] * k[2] - p[2] * k[1] -- not integer
        ----
        So let n be primes p.
        p[i] * k[j] - p[j] * k[i] -- not integer
        ----
        Further example for 3.
        1 * k[1] - 2 * k[0] -- not integer
        1 * k[2] - 3 * k[0] -- not integer
        2 * k[2] - 3 * k[1] -- not integer
        ----
        Symmetry cancels half of combinations, so assume j > i.
        Require k[0] = 0; solve k[j] for k[i].
        ----
        More example for 3.
        k[1] -- not integer
        k[2] -- not integer
        2 * k[2] - 3 * k[1] -- not integer
        Choose thirds.
        1 / 3 -- not integer
        2 / 3 -- not integer
        4 / 3 - 1 -- not integer
        Works, but only for dims >= 3.
        ----
        Notice that a multiset of 1s is coprime.
        Exploiting this allows solution for dim = 2.
        There is no solution for dim = 1.
        */
        double const s = sin(ljfreq(idim) * t + ljphase(idim)) / 2;

        size_div_t const z = size_div(m, n);
        napkin.R[ipoly].r[ibead].d[idim] = ((double) z.rem + (1 + s) / 2) * v;
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
