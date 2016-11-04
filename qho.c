#include "floating.h"
#include "report.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
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

bool qho(void) {
  struct ensemble ens = {
  };
  struct napkin nk = {
  };

  /* p. 44 */

  gsl_rng_env_setup();

  gsl_rng* const rng = gsl_rng_alloc(gsl_rng_mt19937);

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
