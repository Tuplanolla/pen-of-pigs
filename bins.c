#include "bins.h"
#include <math.h>
#include <stdlib.h>

struct bins {
  size_t* m;
  size_t N;
};

void bins_free(struct bins* const bins) {
  free(bins);
}

void bins_forget(struct bins* const bins) {
  bins->N = 0;
}

struct bins* bins_alloc(void) {
  struct bins* const bins = malloc(sizeof *bins);
  if (bins == NULL)
    return NULL;

  bins_forget(bins);

  return bins;
}

void bins_accum(struct bins* const bins, double const x) {}
