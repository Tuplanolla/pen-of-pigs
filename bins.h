#ifndef BINS_H
#define BINS_H

#include "exts.h"
#include <stddef.h>

struct bins;

// The call `bins_free(bins)` releases the memory
// used by the histogram `bins`.
void bins_free(struct bins*);

// The call `bins_forget(bins)` clears the histogram `bins`.
__attribute__ ((__nonnull__))
void bins_forget(struct bins*);

// The call `bins = bins_alloc()` allocates memory
// for the new histogram `bins`.
// When `bins` is no longer used, `bins_free` should be called on it.
__attribute__ ((__malloc__))
struct bins* bins_alloc(void);

// The call `bins_accum(bins, x)` adds `x` to the histogram `bins`.
__attribute__ ((__nonnull__))
void bins_accum(struct bins*, double);

#endif
