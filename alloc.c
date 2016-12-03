#include "alloc.h"
#include <gsl/gsl_rng.h>

void* alloc_err(__attribute__ ((__unused__)) size_t const n) {
  return NULL;
}

static gsl_rng* rng = NULL;
static double p = 0.5;

void alloc_ran_rel(double const q) {
  p = q;
}

void* alloc_ran(size_t const n) {
  if (rng == NULL)
    rng = gsl_rng_alloc(gsl_rng_env_setup());

  if (rng == NULL || gsl_rng_uniform(rng) >= p) {
    errno = ENOMEM;

    return NULL;
  } else
    return malloc(n);
}
