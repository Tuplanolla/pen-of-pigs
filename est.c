#include "est.h"
#include <math.h>
#include <stdlib.h>

struct est {
  size_t N;
  double M1;
  double M2;
};

void est_free(struct est* const est) {
  free(est);
}

struct est* est_alloc(void) {
  struct est* const est = malloc(sizeof *est);
  if (est == NULL)
    return NULL;

  est->N = 0;
  est->M1 = 0.0;
  est->M2 = 0.0;

  return est;
}

void est_accum(struct est* const est, double const x) {
  ++est->N;
  double const dx = x - est->M1;
  est->M1 += dx / (double) est->N;
  est->M2 += dx * (x - est->M1);
}

double est_mean(struct est const* const est) {
  switch (est->N) {
    case 0:
      return NAN;
      break;
    case 1:
      return est->M1;
      break;
    default:
      return est->M1;
  }
}

double est_var(struct est const* const est) {
  switch (est->N) {
    case 0:
      return NAN;
      break;
    case 1:
      return 0.0;
      break;
    default:
      return est->M2 / (double) (est->N - 1);
  }
}

double est_sd(struct est const* const est) {
  switch (est->N) {
    case 0:
      return NAN;
      break;
    case 1:
      return 0.0;
      break;
    default:
      return sqrt(est_var(est));
  }
}

double est_sem(struct est const* const est) {
  switch (est->N) {
    case 0:
      return NAN;
      break;
    case 1:
      return 0.0;
      break;
    default:
      return est_sd(est) / sqrt((double) est->N);
  }
}
