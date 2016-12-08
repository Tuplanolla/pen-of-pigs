#include "exts.h"
#include "lims.h"
#include "sim.h"
#include <errno.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

__attribute__ ((__nonnull__))
static bool qho_parse_n(size_t* const n, size_t const min, size_t const max) {
  char* endptr;
  errno = 0;
  unsigned long long int const k = strtoull(optarg, &endptr, 10);
  if (errno != 0 || endptr == optarg || *endptr != '\0') {
    (void) fprintf(stderr, "Argument '%s' is not an integer.\n", optarg);
    return false;
  } else if (k > SIZE_MAX || (size_t) k < min || (size_t) k > max) {
    (void) fprintf(stderr, "Argument '%s' is out of range.\n", optarg);
    return false;
  } else
    *n = (size_t) k;

  return true;
}

__attribute__ ((__nonnull__))
static bool qho_parse_x(double* const x, double const min, double const max) {
  char* endptr;
  errno = 0;
  double const y = strtod(optarg, &endptr);
  if (errno != 0 || endptr == optarg || *endptr != '\0') {
    (void) fprintf(stderr, "Argument '%s' is not an integer.\n", optarg);
    return false;
  } else if (y < min || y > max) {
    (void) fprintf(stderr, "Argument '%s' is out of range.\n", optarg);
    return false;
  } else
    *x = y;

  return true;
}

__attribute__ ((__nonnull__))
static bool qho_parse(int const n, char** const x,
    size_t* const ndim, size_t* const npoly,
    size_t* const nbead, size_t* const nsubdiv,
    size_t* const nthrm, size_t* const nprod,
    size_t* const nthrmrec, size_t* const nprodrec,
    double* const T) {
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
    int const opt = getopt(n, x, "d:N:M:K:t:p:T:P:h:");
    if (opt == -1)
      break;

    switch ((char) opt) {
      case 'd':
        if (!qho_parse_n(ndim, 1, DIM_MAX))
          return false;
        break;
      case 'N':
        if (!qho_parse_n(npoly, 1, POLY_MAX))
          return false;
        break;
      case 'M':
        if (!qho_parse_n(nbead, 1, BEAD_MAX))
          return false;
        break;
      case 'K':
        if (!qho_parse_n(nsubdiv, 1, SUBDIV_MAX))
          return false;
        break;
      case 't':
        if (!qho_parse_n(nthrm, 0, STEP_MAX))
          return false;
        break;
      case 'p':
        if (!qho_parse_n(nprod, 0, STEP_MAX))
          return false;
        break;
      case 'T':
        if (!qho_parse_n(nthrmrec, 0, STEP_MAX))
          return false;
        break;
      case 'P':
        if (!qho_parse_n(nprodrec, 0, STEP_MAX))
          return false;
        break;
      case 'h':
        if (!qho_parse_x(T, 0.0, INFINITY))
          return false;
        break;
      default:
        (void) fprintf(stderr, "Argument '%s' is invalid.\n", optarg);
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
        &nthrm, &nprod, &nthrmrec, &nprodrec, &T))
    return EXIT_FAILURE;

  double const beta = 1.0 / T;

  double const q = omega / 2.0;
  double const E = (double) ndim * q / tanh(q * beta);
  (void) printf("Expected for QHO: E = %f (T = %f)\n", E, T);

  return sim_run(ndim, npoly, nbead, nsubdiv, nthrm, nprod, nthrmrec, nprodrec,
      false, 10.0, beta, pot_lj612, pot_zero, potext_harm) ?
    EXIT_SUCCESS : EXIT_FAILURE;
}
