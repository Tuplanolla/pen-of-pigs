#include "hist.h"
#include "fp.h"
#include "size.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

struct hist {
  size_t* m;
  size_t* tmp;
  size_t ndim;
  size_t nlin;
  size_t ncub;
  double a;
  double b;
  size_t nsum;
};

void hist_free(struct hist* const hist) {
  if (hist != NULL) {
    free(hist->tmp);

    free(hist->m);
  }

  free(hist);
}

void hist_forget(struct hist* const hist) {
  for (size_t i = 0; i < hist->ncub; ++i)
    hist->m[i] = 0;

  hist->nsum = 0;
}

struct hist* hist_alloc(size_t const ndim, size_t const nlin,
    double const a, double const b) {
  bool p = true;

  struct hist* const hist = malloc(sizeof *hist);
  if (hist == NULL)
    p = false;
  else {
    hist->ndim = ndim;
    hist->nlin = nlin;
    hist->ncub = size_pow(nlin, ndim);
    hist->a = a;
    hist->b = b;

    hist->m = malloc(hist->ncub * sizeof *hist->m);
    if (hist->m == NULL)
      p = false;
    else
      hist_forget(hist);

    hist->tmp = malloc(ndim * sizeof *hist->tmp);
    if (hist->tmp == NULL)
      p = false;
  }

  if (p)
    return hist;
  else {
    hist_free(hist);

    return NULL;
  }
}

size_t hist_bin(struct hist const* const hist, double const* const x) {
  for (size_t i = 0; i < hist->ndim; ++i)
    if (x[i] >= hist->a && x[i] < hist->b)
      hist->tmp[i] = size_uclamp((size_t) floor(fp_lerp(x[i],
              hist->a, hist->b, 0.0, (double) hist->nlin)), hist->nlin);
    else
      return SIZE_MAX;

  return size_unhc(hist->nlin, hist->tmp, hist->ndim);
}

__attribute__ ((__nonnull__))
static void unbin(struct hist const* const hist, double* const x,
    size_t const j, double const h) {
  size_hc(j, hist->nlin, hist->tmp, hist->ndim);

  for (size_t i = 0; i < hist->ndim; ++i)
    x[i] = fp_lerp((double) hist->tmp[i] + h,
        0.0, (double) hist->nlin, hist->a, hist->b);
}

void hist_unbin(struct hist const* const hist, double* const x,
    size_t const j) {
  unbin(hist, x, j, 0.5);
}

void hist_funbin(struct hist const* const hist, double* const x,
    size_t const j) {
  unbin(hist, x, j, 0.0);
}

void hist_cunbin(struct hist const* const hist, double* const x,
    size_t const j) {
  unbin(hist, x, j, 1.0);
}

bool hist_accum(struct hist* const hist, double const* const x) {
  size_t const i = hist_bin(hist, x);

  if (i == SIZE_MAX)
    return false;
  else {
    ++hist->m[i];

    ++hist->nsum;

    return true;
  }
}

size_t hist_ndim(struct hist const* const hist) {
  return hist->ndim;
}

size_t hist_nsubdiv(struct hist const* const hist) {
  return hist->nlin;
}

size_t hist_nbin(struct hist const* const hist) {
  return hist->ncub;
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
  return hist->nsum;
}

double hist_normhits(struct hist const* const hist, size_t const i) {
  return hist->nsum == 0 ? 0.0 : (double) hist->m[i] / (double) hist->nsum;
}
