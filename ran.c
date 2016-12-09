#include "exts.h"
#include "ran.h"
#include <gsl/gsl_rng.h>
#include <stddef.h>
#include <stdio.h>
#include <time.h>

size_t ran_index(gsl_rng* const rng, size_t const n) {
  return (size_t) gsl_rng_uniform_int(rng, n);
}

double ran_uopen(gsl_rng* const rng, double const x) {
  return gsl_rng_uniform_pos(rng) * x;
}

double ran_open(gsl_rng* const rng, double const x) {
  return (2.0 * gsl_rng_uniform_pos(rng) - 1.0) * x;
}

__attribute__ ((__const__, __pure__))
static unsigned long int filog16(unsigned long int n) {
  size_t m = 0;

  while (n >= 16) {
    n /= 16;
    ++m;
  }

  return m;
}

bool ran_dateid(gsl_rng* const rng, char const* const pre,
    char* const buf, size_t const n) {
  time_t t;
  if (time(&t) == (time_t) -1)
    return false;

  struct tm tm;
  if (localtime_r(&t, &tm) == NULL)
    return false;

  int const k = snprintf(buf, n,
      "%s-%04d-%02d-%02d-%02d-%02d-%02d-%0*lx", pre,
      tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
      tm.tm_hour, tm.tm_min, tm.tm_sec,
      1 + (int) filog16(gsl_rng_max(rng)), gsl_rng_get(rng));
  if (k < 0 || (size_t) k >= n)
    return false;

  return true;
}
