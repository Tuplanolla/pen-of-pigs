#ifndef HIST_H
#define HIST_H

#include "exts.h"
#include <stddef.h>

// This structure holds hit counts and cached values.
struct hist;

// The call `hist_free(hist)` releases the memory
// used by the histogram `hist`.
void hist_free(struct hist*);

// The call `hist_forget(hist)` clears the histogram `hist`.
__attribute__ ((__nonnull__))
void hist_forget(struct hist*);

// The statement `hist = hist_alloc(ndim, nlin, a, b)` allocates memory
// for managing the histogram `hist` of `ndim` dimensions,
// `nlin` bins per dimension and a range from `a` to `b`.
// When `hist` is no longer used, `hist_free` should be called on it.
__attribute__ ((__malloc__))
struct hist* hist_alloc(size_t, size_t, double, double);

// The call `hist_bin(hist, x)` returns
// the index of the bin for `x` in the histogram `hist`
// if `x` is in the valid range.
// Otherwise `SIZE_MAX` is returned.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t hist_bin(struct hist const*, double const*);

// The call `hist_unbin(hist, x, i)` writes
// the center of the bin `i` from the histogram `hist` into `x`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
void hist_unbin(struct hist const*, double*, size_t);

// The call `hist_funbin(hist, x, i)` writes
// the minimum of the bin `i` from the histogram `hist` into `x`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
void hist_funbin(struct hist const*, double*, size_t);

// The call `hist_cunbin(hist, x, i)` writes
// the maximum of the bin `i` from the histogram `hist` into `x`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
void hist_cunbin(struct hist const*, double*, size_t);

// The call `hist_accum(hist, x)` adds `x` to the histogram `hist` and
// returns `true` if `x` is in the valid range.
// Otherwise `false` is returned.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
bool hist_accum(struct hist*, double const*);

// The call `hist_ndim(hist)` returns
// the number of dimensions in the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t hist_ndim(struct hist const*);

// The call `hist_nsubdiv(hist)` returns
// the number of bins per dimension in the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t hist_nsubdiv(struct hist const*);

// The call `hist_nbin(hist)` returns
// the number of bins in the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t hist_nbin(struct hist const*);

// The call `hist_min(hist)` returns
// the lower limit of the valid range of the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double hist_min(struct hist const*);

// The call `hist_max(hist)` returns
// the higher limit of the valid range of the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double hist_max(struct hist const*);

// The call `hist_length(hist)` returns
// the length of one dimension of the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double hist_length(struct hist const*);

// The call `hist_volume(hist)` returns
// the product of the length of all the dimensions of the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double hist_volume(struct hist const*);

// The call `hist_hits(hist, i)` returns
// the hit count for bin `i` in the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t hist_hits(struct hist const*, size_t);

// The call `hist_sumhits(hist)` returns
// the total hit count in the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t hist_sumhits(struct hist const*);

// The call `hist_normhits(hist, i)` returns
// the normalized hit count for bin `i` in the histogram `hist`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double hist_normhits(struct hist const*, size_t);

#endif
