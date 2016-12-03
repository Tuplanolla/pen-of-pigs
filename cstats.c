#include "cstats.h"
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <stdlib.h>

struct cstats {
  double* x;
  size_t M;
  size_t N;
  double M1;
  double M2;
};

void cstats_free(struct cstats* const cstats) {
  if (cstats != NULL)
    free(cstats->x);

  free(cstats);
}

void cstats_forget(struct cstats* const cstats) {
  cstats->N = 0;
  cstats->M1 = 0.0;
  cstats->M2 = 0.0;
}

struct cstats* cstats_alloc(size_t const M) {
  bool e = false;

  struct cstats* const cstats = malloc(sizeof *cstats);
  if (cstats == NULL)
    e = true;
  else {
    cstats->x = malloc(M * sizeof *cstats->x);
    if (cstats == NULL)
      e = true;
    else
      cstats->M = M;
  }

  cstats_forget(cstats);

  if (e) {
    cstats_free(cstats);

    return NULL;
  } else
    return cstats;
}

bool cstats_accum(struct cstats* const cstats, double const x) {
  if (cstats->N < cstats->M) {
    cstats->x[cstats->N] = x;
    ++cstats->N;
    double const dx = x - cstats->M1;
    cstats->M1 += dx / (double) cstats->N;
    cstats->M2 += dx * (x - cstats->M1);

    return true;
  } else
    return false;
}

double cstats_mean(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    default:
      return cstats->M1;
  }
}

double cstats_var(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return cstats->M2 / (double) (cstats->N - 1);
  }
}

double cstats_sd(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(cstats_var(cstats));
  }
}

double cstats_autocorr(struct cstats const* const cstats, size_t const k) {
  double S = 0.0;

  switch (cstats->N) {
    case 0:
    case 1:
      return NAN;
    default:
      for (size_t i = k; i < cstats->N; ++i)
        S += (cstats->x[i - k] - cstats->M1) * (cstats->x[i] - cstats->M1);

      return ((double) (cstats->N - 1) / (double) (cstats->N - k)) *
        (S / cstats->M2);
  }
}

double cstats_corrtime(struct cstats const* const cstats) {
  double S = 0.0;

  switch (cstats->N) {
    case 0:
    case 1:
      return NAN;
    default:
      for (size_t k = 1; k < cstats->N; ++k) {
        double const R = cstats_autocorr(cstats, k);
        if (R < 0.0)
          break;

        S += R;
      }

      return 1.0 + 2.0 * S;
  }
}

double cstats_sem(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(cstats_var(cstats) / (double) cstats->N);
  }
}

double cstats_corrsem(struct cstats const* const cstats) {
  switch (cstats->N) {
    case 0:
      return NAN;
    case 1:
      return 0.0;
    default:
      return sqrt(cstats_corrtime(cstats) *
          cstats_var(cstats) / (double) cstats->N);
  }
}
