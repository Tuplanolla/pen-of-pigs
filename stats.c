#include "stats.h"
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <stdlib.h>

struct stats {
  double* x;
  size_t M;
  size_t N;
  double M1;
  double M2;
};

void stats_free(struct stats* const stats) {
  if (stats != NULL)
    free(stats->x);

  free(stats);
}

void stats_forget(struct stats* const stats) {
  stats->N = 0;
  stats->M1 = 0.0;
  stats->M2 = 0.0;
}

struct stats* stats_alloc(size_t const M) {
  bool e = false;

  struct stats* const stats = malloc(sizeof *stats);
  if (stats == NULL)
    e = true;
  else {
    stats->x = malloc(M * sizeof *stats->x);
    if (stats == NULL)
      e = true;
    else
      stats->M = M;
  }

  stats_forget(stats);

  if (e) {
    stats_free(stats);

    return NULL;
  } else
    return stats;
}

bool stats_accum(struct stats* const stats, double const x) {
  if (stats->N < stats->M) {
    stats->x[stats->N] = x;
    ++stats->N;
    double const dx = x - stats->M1;
    stats->M1 += dx / (double) stats->N;
    stats->M2 += dx * (x - stats->M1);

    return true;
  } else
    return false;
}

double stats_mean(struct stats const* const stats) {
  switch (stats->N) {
    case 0:
      return NAN;
    default:
      return stats->M1;
  }
}

double stats_var(struct stats const* const stats) {
  switch (stats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return stats->M2 / (double) (stats->N - 1);
  }
}

double stats_sd(struct stats const* const stats) {
  switch (stats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_var(stats));
  }
}

double stats_sem(struct stats const* const stats) {
  switch (stats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_var(stats) / (double) stats->N);
  }
}

double stats_autocorr(struct stats const* const stats, size_t const k) {
  double S = 0.0;

  switch (stats->N) {
    case 0:
    case 1:
      return NAN;
    default:
      for (size_t i = k; i < stats->N; ++i)
        S += (stats->x[i - k] - stats->M1) * (stats->x[i] - stats->M1);

      return ((double) (stats->N - 1) / (double) (stats->N - k)) *
        (S / stats->M2);
  }
}

double stats_corrtime(struct stats const* const stats) {
  double S = 0.0;

  switch (stats->N) {
    case 0:
    case 1:
      return NAN;
    default:
      for (size_t k = 1; k < stats->N; ++k) {
        double const R = stats_autocorr(stats, k);
        if (R < 0.0)
          break;

        S += R;
      }

      return 1.0 + 2.0 * S;
  }
}

double stats_corrsem(struct stats const* const stats) {
  switch (stats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(stats_corrtime(stats) *
          stats_var(stats) / (double) stats->N);
  }
}
