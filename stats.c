#include "stats.h"
#include <math.h>
#include <stdlib.h>

struct stats {
  size_t N;
  double M1;
  double M2;
};

void stats_free(struct stats* const stats) {
  free(stats);
}

void stats_forget(struct stats* const stats) {
  stats->N = 0;
  stats->M1 = 0.0;
  stats->M2 = 0.0;
}

struct stats* stats_alloc(void) {
  struct stats* const stats = malloc(sizeof *stats);
  if (stats == NULL)
    return NULL;

  stats_forget(stats);

  return stats;
}

void stats_accum(struct stats* const stats, double const x) {
  ++stats->N;
  double const dx = x - stats->M1;
  stats->M1 += dx / (double) stats->N;
  stats->M2 += dx * (x - stats->M1);
}

double stats_mean(struct stats const* const stats) {
  switch (stats->N) {
    case 0:
      return NAN;
    case 1:
      return stats->M1;
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
