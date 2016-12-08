#include "exts.h"
#include "lims.h"
#include "opt.h"
#include "sim.h"
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

__attribute__ ((__nonnull__))
static bool he4_parse(int const n, char* const* const x,
    size_t* const ndim, size_t* const npoly,
    size_t* const nbead, size_t* const nsubdiv,
    size_t* const nthrm, size_t* const nprod,
    size_t* const nthrmrec, size_t* const nprodrec) {
  char const* shortstr = "d:N:M:K:h:p:H:P:";
  char const* const longstr[] = {
    "ndim", "npoly", "nbead", "nsubdiv",
    "nthrm", "nprod", "nthrmrec", "nprodrec"
  };

  *ndim = 1;
  *npoly = 1;
  *nbead = 1;
  *nsubdiv = 1;
  *nthrm = 0;
  *nprod = 0;
  *nthrmrec = 0;
  *nprodrec = 0;

  for ever {
    int const i = opt_parse(n, x, shortstr, longstr);
    if (i == -1)
      break;

    switch ((char) i) {
      case 'd':
        if (!opt_parse_size(ndim, 1, DIM_MAX))
          return false;
        break;
      case 'N':
        if (!opt_parse_size(npoly, 1, POLY_MAX))
          return false;
        break;
      case 'M':
        if (!opt_parse_size(nbead, 1, BEAD_MAX))
          return false;
        break;
      case 'K':
        if (!opt_parse_size(nsubdiv, 1, SUBDIV_MAX))
          return false;
        break;
      case 'h':
        if (!opt_parse_size(nthrm, 0, STEP_MAX))
          return false;
        break;
      case 'p':
        if (!opt_parse_size(nprod, 0, STEP_MAX))
          return false;
        break;
      case 'H':
        if (!opt_parse_size(nthrmrec, 0, STEP_MAX))
          return false;
        break;
      case 'P':
        if (!opt_parse_size(nprodrec, 0, STEP_MAX))
          return false;
        break;
      default:
        return false;
    }
  }

  return true;
}

static double const epsilon = 4.0;
static double const sigma = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double pot_lj612(struct ensem const* const ensem,
    struct bead const* const r0, struct bead const* const r1) {
  double const d = bead_dist2(ensem, r0, r1);

  if (d == 0.0)
    return INFINITY;
  else {
    double const sigmar2 = gsl_pow_2(sigma) / d;

    return 4.0 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
  }
}

int main(int const n, char** const x) {
  size_t ndim;
  size_t npoly;
  size_t nbead;
  size_t nsubdiv;
  size_t nthrm;
  size_t nprod;
  size_t nthrmrec;
  size_t nprodrec;

  if (!he4_parse(n, x, &ndim, &npoly, &nbead, &nsubdiv,
        &nthrm, &nprod, &nthrmrec, &nprodrec)) {
    (void) fprintf(stderr, "Failed to parse argument list.\n");

    return EXIT_FAILURE;
  }

  if (!sim_run(ndim, npoly, nbead, nsubdiv, nthrm, nprod, nthrmrec, nprodrec,
      true, 10.0, 1e+3, pot_lj612, pot_zero, potext_zero)) {
    (void) fprintf(stderr, "Failed to run simulation.\n");

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
