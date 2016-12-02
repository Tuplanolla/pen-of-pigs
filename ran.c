#include "ran.h"
#include <gsl/gsl_rng.h>
#include <stdlib.h>

__attribute__ ((__nonnull__))
size_t ran_index(gsl_rng* const rng, size_t const n) {
  return (size_t) gsl_rng_uniform_int(rng, n);
}

__attribute__ ((__nonnull__))
double ran_uopen(gsl_rng* const rng, double const x) {
  return x * gsl_rng_uniform_pos(rng);
}

__attribute__ ((__nonnull__))
double ran_open(gsl_rng* const rng, double const x) {
  return x * (2.0 * gsl_rng_uniform_pos(rng) - 1);
}
