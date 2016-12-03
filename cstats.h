#ifndef CSTATS_H
#define CSTATS_H

#include "exts.h"
#include <stddef.h>

// This structure holds a cache of values for calculating their statistics.
struct cstats;

// The call `cstats_free(cstats)` releases the memory
// used by the cache `cstats`.
void cstats_free(struct cstats*);

// The call `cstats_forget(cstats)` clears the cache `cstats`.
__attribute__ ((__nonnull__))
void cstats_forget(struct cstats*);

// The call `cstats = cstats_alloc()` allocates memory
// for the new cache `cstats`.
// When `cstats` is no longer used, `cstats_free` should be called on it.
__attribute__ ((__malloc__))
struct cstats* cstats_alloc(size_t);

// The call `cstats_accum(cstats, x)` adds `x` to the cache `cstats`.
__attribute__ ((__nonnull__))
bool cstats_accum(struct cstats*, double);

// The call `cstats_corrtime(cstats)` returns
// the estimated correlation time of the cache `cstats`.
__attribute__ ((__nonnull__, __pure__))
double cstats_corrtime(struct cstats const*);

// The call `cstats_mean(cstats)` returns
// the arithmetic mean of the cache `cstats`.
__attribute__ ((__nonnull__, __pure__))
double cstats_mean(struct cstats const*);

// The call `cstats_var(cstats)` returns the variance of the cache `cstats`.
__attribute__ ((__nonnull__, __pure__))
double cstats_var(struct cstats const*);

// The call `cstats_sd(cstats)` returns
// the standard deviation of the cache `cstats`.
__attribute__ ((__nonnull__, __pure__))
double cstats_sd(struct cstats const*);

// The call `cstats_sem(cstats)` returns
// the standard error of the mean of the cache `cstats`.
__attribute__ ((__nonnull__, __pure__))
double cstats_sem(struct cstats const*);

#endif
