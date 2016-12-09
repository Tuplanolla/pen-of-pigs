#include "stats.h"
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

struct stats {
  double* x;
  size_t N;
  size_t n;
  double M1;
  double M2;
};

void stats_free(struct stats* const stats) {
  if (stats != NULL)
    free(stats->x);

  free(stats);
}

void stats_forget(struct stats* const stats) {
  stats->n = 0;
  stats->M1 = 0.0;
  stats->M2 = 0.0;
}

struct stats* stats_alloc(size_t const N) {
  bool p = true;

  struct stats* const stats = malloc(sizeof *stats);
  if (stats == NULL)
    p = false;
  else if (N >= SIZE_MAX / sizeof *stats->x)
    p = false;
  else {
    stats->x = malloc(N * sizeof *stats->x);
    if (stats->x == NULL)
      p = false;
    else {
      stats->N = N;

      stats_forget(stats);
    }
  }

  if (p)
    return stats;
  else {
    stats_free(stats);

    return NULL;
  }
}

bool stats_accum(struct stats* const stats, double const x) {
  if (stats->n < stats->N) {
    stats->x[stats->n] = x;
    ++stats->n;
    double const dx = x - stats->M1;
    stats->M1 += dx / (double) stats->n;
    stats->M2 += dx * (x - stats->M1);

    return true;
  } else
    return false;
}

size_t stats_n(struct stats const* const stats) {
  return stats->n;
}

double stats_mean(struct stats const* const stats) {
  switch (stats->n) {
    case 0:
      return NAN;
    default:
      return stats->M1;
  }
}

double stats_var(struct stats const* const stats) {
  switch (stats->n) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return stats->M2 / (double) (stats->n - 1);
  }
}

double stats_sd(struct stats const* const stats) {
  switch (stats->n) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_var(stats));
  }
}

double stats_sem(struct stats const* const stats) {
  switch (stats->n) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_var(stats) / (double) stats->n);
  }
}

double stats_autocorr(struct stats const* const stats, size_t const k) {
  double S = 0.0;

  switch (stats->n) {
    case 0:
    case 1:
      return NAN;
    default:
      for (size_t i = k; i < stats->n; ++i)
        S += (stats->x[i - k] - stats->M1) * (stats->x[i] - stats->M1);

      return ((double) (stats->n - 1) / (double) (stats->n - k)) *
        (S / stats->M2);
  }
}

double stats_corrtime(struct stats const* const stats) {
  double S = 0.0;

  switch (stats->n) {
    case 0:
    case 1:
      return NAN;
    default:
      for (size_t k = 1; k < stats->n; ++k) {
        double const A = stats_autocorr(stats, k);
        if (A < 0.0)
          break;

        S += A;
      }

      return 1.0 + 2.0 * S;
  }
}

double stats_corrsem(struct stats const* const stats) {
  switch (stats->n) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_corrtime(stats) *
          stats_var(stats) / (double) stats->n);
  }
}
