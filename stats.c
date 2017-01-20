#include "size.h"
#include "stats.h"
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

struct stats {
  double* x;
  size_t ncap;
  size_t nmemb;
  double m1;
  double m2;
};

void stats_free(struct stats* const stats) {
  if (stats != NULL)
    free(stats->x);

  free(stats);
}

void stats_forget(struct stats* const stats) {
  stats->nmemb = 0;
  stats->m1 = 0.0;
  stats->m2 = 0.0;
}

struct stats* stats_alloc(size_t const ncap, bool const cache) {
  bool p = true;

  struct stats* const stats = malloc(sizeof *stats);
  if (stats == NULL)
    p = false;
  else {
    stats->ncap = ncap;

    if (cache) {
      stats->x = malloc(ncap * sizeof *stats->x);
      if (stats->x == NULL)
        p = false;
    } else
      stats->x = NULL;

    stats_forget(stats);
  }

  if (p)
    return stats;
  else {
    stats_free(stats);

    return NULL;
  }
}

bool stats_accum(struct stats* const stats, double const x) {
  if (stats->nmemb < stats->ncap) {
    if (stats->x != NULL)
      stats->x[stats->nmemb] = x;

    ++stats->nmemb;

    double const dx = x - stats->m1;
    stats->m1 += dx / (double) stats->nmemb;
    stats->m2 += dx * (x - stats->m1);

    return true;
  } else
    return false;
}

bool stats_cache(struct stats const* const stats) {
  return stats->x != NULL;
}

size_t stats_n(struct stats const* const stats) {
  return stats->nmemb;
}

double stats_mean(struct stats const* const stats) {
  switch (stats->nmemb) {
    case 0:
      return (double) NAN;
    default:
      return stats->m1;
  }
}

double stats_var(struct stats const* const stats) {
  switch (stats->nmemb) {
    case 0:
      return (double) NAN;
    case 1:
      return 0.0;
    default:
      return stats->m2 / (double) (stats->nmemb - 1);
  }
}

double stats_sd(struct stats const* const stats) {
  switch (stats->nmemb) {
    case 0:
      return (double) NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_var(stats));
  }
}

double stats_sem(struct stats const* const stats) {
  switch (stats->nmemb) {
    case 0:
      return (double) NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_var(stats) / (double) stats->nmemb);
  }
}

double stats_autocorr(struct stats const* const stats, size_t const k) {
  if (stats->x == NULL)
    return (double) NAN;

  double s = 0.0;

  switch (stats->nmemb) {
    case 0:
    case 1:
      return (double) NAN;
    default:
      for (size_t i = k; i < stats->nmemb; ++i)
        s += (stats->x[i - k] - stats->m1) * (stats->x[i] - stats->m1);

      return ((double) (stats->nmemb - 1) / (double) (stats->nmemb - k)) *
        (s / stats->m2);
  }
}

double stats_corrtime_lag(struct stats const* const stats, size_t n) {
  if (stats->x == NULL)
    return (double) NAN;

  double s = 0.0;

  switch (stats->nmemb) {
    case 0:
    case 1:
      return (double) NAN;
    default:
      n = size_min(n, stats->nmemb);

      for (size_t k = 1; k < n; ++k) {
        double const a = stats_autocorr(stats, k);
        if (a >= 0.0)
          s += a;
        else
          break;
      }

      return 1.0 + 2.0 * s;
  }
}

double stats_corrtime(struct stats const* const stats) {
  return stats_corrtime_lag(stats, SIZE_MAX);
}

double stats_corrsem(struct stats const* const stats) {
  if (stats->x == NULL)
    return (double) NAN;

  switch (stats->nmemb) {
    case 0:
      return (double) NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_corrtime(stats) *
          stats_var(stats) / (double) stats->nmemb);
  }
}
