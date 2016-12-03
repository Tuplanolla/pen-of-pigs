#include "alloc.h"
#include "cstats.h"
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <stdlib.h>

static double autocorr(double const* const x, size_t const n,
    double const mean, double const var, size_t const k) {
  double s = 0.0;

  for (size_t i = k; i < n; ++i)
    s += (x[i - k] - mean) * (x[i] - mean);

  return s / ((double) (n - k) * var);
}

static double corrtime(double const* const x, size_t const n,
    double const mean, double const var) {
  double s = 0.0;

  for (size_t k = 1; k < n; ++k) {
    double const C = autocorr(x, n, mean, var, k);
    if (C < 0)
      break;

    s += C;
  }

  return 1.0 + 2.0 * s;
}

struct cstats {
  size_t M;
  double* x;
  size_t N;
};

void cstats_free(struct cstats* const cstats) {
  if (cstats != NULL)
    free(cstats->x);

  free(cstats);
}

void cstats_forget(struct cstats* const cstats) {
  cstats->N = 0;
}

struct cstats* cstats_alloc(size_t const M) {
  alloc_proc alloc = malloc;

  struct cstats* const cstats = alloc(sizeof *cstats);
  if (cstats == NULL)
    alloc = alloc_err;
  else {
    cstats->M = M;

    cstats->x = malloc(cstats->M * sizeof *cstats->x);
    if (cstats == NULL)
      alloc = alloc_err;
  }

  cstats_forget(cstats);

  if (alloc == alloc_err) {
    cstats_free(cstats);

    return NULL;
  } else
    return cstats;
}

bool cstats_accum(struct cstats* const cstats, double const x) {
  if (cstats->N < cstats->M) {
    cstats->x[cstats->N] = x;
    ++cstats->N;

    return true;
  } else
    return false;
}

double cstats_corrtime(struct cstats const* const cstats) {
  return corrtime(cstats->x, cstats->N, cstats_mean(cstats), cstats_var(cstats));
}

double cstats_mean(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    case 1:
      return cstats->x[0];
    default:
      return gsl_stats_mean(cstats->x, 1, cstats->N);
  }
}

double cstats_var(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return gsl_stats_variance(cstats->x, 1, cstats->N);
  }
}

double cstats_sd(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return gsl_stats_sd(cstats->x, 1, cstats->N);
  }
}

double cstats_sem(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default: {
        double const var = cstats_var(cstats);
        return sqrt(corrtime(cstats->x, cstats->N, cstats_mean(cstats), var) *
            var / (double) cstats->N);
      }
  }
}
