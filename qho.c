#include "exts.h"
#include "lims.h"
#include "opt.h"
#include "sim.h"
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

__attribute__ ((__nonnull__))
static bool qho_parse(int const n, char* const* const x,
    size_t* const ndim, size_t* const npoly,
    size_t* const nbead, size_t* const nsubdiv,
    size_t* const nthrm, size_t* const nprod,
    size_t* const nthrmrec, size_t* const nprodrec,
    double* const T) {
  char const* shortstr = "d:N:M:K:h:p:H:P:T:";
  char const* const longstr[] = {
    "ndim", "npoly", "nbead", "nsubdiv",
    "nthrm", "nprod", "nthrmrec", "nprodrec",
    "temp"
  };

  *ndim = 1;
  *npoly = 1;
  *nbead = 1;
  *nsubdiv = 1;
  *nthrm = 0;
  *nprod = 0;
  *nthrmrec = 0;
  *nprodrec = 0;
  *T = 1.0;

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
      case 'T':
        if (!opt_parse_fp(T, 0.0, INFINITY))
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
  double const sigmar2 = gsl_pow_2(sigma) / bead_dist2(ensem, r0, r1);

  return 4.0 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

static double const omega = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double potext_harm(struct ensem const* const ensem,
    struct bead const* const r) {
  return gsl_pow_2(omega) * bead_norm2(ensem, r) / 2.0;
}

__attribute__ ((__nonnull__))
int main(int const n, char** const x) {
  size_t ndim;
  size_t npoly;
  size_t nbead;
  size_t nsubdiv;
  size_t nthrm;
  size_t nprod;
  size_t nthrmrec;
  size_t nprodrec;
  double T;

  if (!qho_parse(n, x, &ndim, &npoly, &nbead, &nsubdiv,
        &nthrm, &nprod, &nthrmrec, &nprodrec, &T)) {
    (void) fprintf(stderr, "Failed to parse argument list.\n");

    return EXIT_FAILURE;
  }

  double const beta = 1.0 / T;

  double const q = omega / 2.0;
  double const E = (double) ndim * q / tanh(q * beta);
  (void) printf("Expected for QHO: E = %f (T = %f)\n", E, T);

  if (!sim_run(ndim, npoly, nbead, nsubdiv, nthrm, nprod, nthrmrec, nprodrec,
      false, 10.0, beta, pot_zero, pot_zero, potext_harm)) {
    (void) fprintf(stderr, "Failed to run simulation.\n");

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
