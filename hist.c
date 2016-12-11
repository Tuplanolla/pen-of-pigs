#include "hist.h"
#include "fp.h"
#include "size.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

struct hist {
  size_t* m;
  size_t* t;
  size_t d;
  size_t N;
  size_t M;
  double a;
  double b;
  size_t n;
};

void hist_free(struct hist* const hist) {
  if (hist != NULL) {
    free(hist->t);

    free(hist->m);
  }

  free(hist);
}

void hist_forget(struct hist* const hist) {
  for (size_t i = 0; i < hist->M; ++i)
    hist->m[i] = 0;

  hist->n = 0;
}

struct hist* hist_alloc(size_t const d, size_t const N,
    double const a, double const b) {
  bool p = true;

  struct hist* const hist = malloc(sizeof *hist);
  if (hist == NULL)
    p = false;
  else {
    hist->d = d;
    hist->N = N;
    hist->M = size_pow(N, d);
    hist->a = a;
    hist->b = b;

    hist->m = malloc(hist->M * sizeof *hist->m);
    if (hist->m == NULL)
      p = false;
    else
      hist_forget(hist);

    hist->t = malloc(d * sizeof *hist->t);
    if (hist->t == NULL)
      p = false;
  }

  if (p)
    return hist;
  else {
    hist_free(hist);

    return NULL;
  }
}

size_t hist_bindex(struct hist const* const hist, double const* const x) {
  for (size_t i = 0; i < hist->d; ++i)
    if (x[i] >= hist->a && x[i] < hist->b)
      hist->t[i] = size_uclamp((size_t) floor(fp_lerp(x[i],
              hist->a, hist->b, 0.0, (double) hist->N)), hist->N);
    else
      return SIZE_MAX;

  return size_unhc(hist->N, hist->t, hist->d);
}

void hist_unbindex(struct hist const* const hist, double* const x,
    size_t const j) {
  size_hc(j, hist->M, hist->t, hist->d);

  for (size_t i = 0; i < hist->d; ++i)
    x[i] = fp_lerp((double) hist->t[i] + 0.5,
        0.0, (double) hist->N, hist->a, hist->b);
}

bool hist_accum(struct hist* const hist, double const* const x) {
  size_t const i = hist_bindex(hist, x);

  if (i == SIZE_MAX)
    return false;
  else {
    ++hist->m[i];

    ++hist->n;

    return true;
  }
}

size_t hist_ndim(struct hist const* const hist) {
  return hist->d;
}

size_t hist_nsubdiv(struct hist const* const hist) {
  return hist->N;
}

size_t hist_nbin(struct hist const* const hist) {
  return hist->M;
}

double hist_min(struct hist const* const hist) {
  return hist->a;
}

double hist_max(struct hist const* const hist) {
  return hist->b;
}

size_t hist_hits(struct hist const* const hist, size_t const i) {
  return hist->m[i];
}

size_t hist_sumhits(struct hist const* const hist) {
  return hist->n;
}

double hist_normhits(struct hist const* const hist, size_t const i) {
  return hist->n == 0 ? 0.0 : (double) hist->m[i] / (double) hist->n;
}
