#ifndef STATS_H
#define STATS_H

#include "exts.h"
#include <stddef.h>

// This structure holds samples and their statistics.
struct stats;

// The call `stats_free(stats)` releases the memory
// used by the statistics `stats`.
void stats_free(struct stats*);

// The call `stats_forget(stats)` clears the statistics `stats`.
__attribute__ ((__nonnull__))
void stats_forget(struct stats*);

// The statement `stats = stats_alloc(n, p)` allocates memory
// for tracking the statistics `stats` of at most `n` samples,
// caching the samples if `p` is true.
// Otherwise `NULL` is returned.
// When `stats` is no longer used, `stats_free` should be called on it.
__attribute__ ((__malloc__))
struct stats* stats_alloc(size_t, bool);

// The call `stats_accum(stats, x)` adds `x` to the statistics `stats` and
// returns `true` if the maximum number of samples has not been reached.
// Otherwise `false` is returned.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
bool stats_accum(struct stats*, double);

// The call `stats_cache(stats)` returns
// `true` if the statistics `stats` caches its samples.
// Otherwise `false` is returned.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
bool stats_cache(struct stats const*);

// The call `stats_n(stats)` returns
// the number of samples in the statistics `stats`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t stats_n(struct stats const*);

// The call `stats_mean(stats)` returns
// the arithmetic mean of the statistics `stats`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double stats_mean(struct stats const*);

// The call `stats_var(stats)` returns the variance of the statistics `stats`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double stats_var(struct stats const*);

// The call `stats_sd(stats)` returns
// the standard deviation of the statistics `stats`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double stats_sd(struct stats const*);

// The call `stats_sem(stats)` returns
// the standard error of the mean of the statistics `stats`.
// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double stats_sem(struct stats const*);

// The call `stats_autocorr(stats, k)` returns
// the estimated lag `k` autocorrelation of the statistics `stats`
// if `stats` caches its samples.
// Otherwise `NAN` is returned.
// The time complexity is $O(n)$ for $n$ samples.
__attribute__ ((__nonnull__, __pure__))
double stats_autocorr(struct stats const*, size_t);

// The call `stats_corrtime(stats)` returns
// the estimated correlation time of the statistics `stats`
// if `stats` caches its samples.
// Otherwise `NAN` is returned.
// The time complexity is $O(n^2)$ for $n$ samples.
__attribute__ ((__nonnull__, __pure__))
double stats_corrtime(struct stats const*);

// The call `stats_corrsem(stats)` returns
// the correlated standard error of the mean of the statistics `stats`
// if `stats` caches its samples.
// Otherwise `NAN` is returned.
// The time complexity is $O(n^2)$ for $n$ samples.
__attribute__ ((__nonnull__, __pure__))
double stats_corrsem(struct stats const*);

#endif
