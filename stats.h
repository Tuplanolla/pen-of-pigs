#ifndef STATS_H
#define STATS_H

#include "exts.h"
#include <stddef.h>

// This structure holds statistics.
struct stats;

// The call `stats_free(stats)` releases the memory
// used by the statistics `stats`.
void stats_free(struct stats*);

// The call `stats_forget(stats)` clears the statistics `stats`.
__attribute__ ((__nonnull__))
void stats_forget(struct stats*);

// The call `stats = stats_alloc()` allocates memory
// for the new statistics `stats`.
// When `stats` is no longer used, `stats_free` should be called on it.
__attribute__ ((__malloc__))
struct stats* stats_alloc(void);

// The call `stats_accum(stats, x)` adds `x` to the statistics `stats`.
__attribute__ ((__nonnull__))
void stats_accum(struct stats*, double);

// The call `stats_mean(stats)` returns
// the arithmetic mean of the statistics `stats`.
__attribute__ ((__nonnull__, __pure__))
double stats_mean(struct stats const*);

// The call `stats_var(stats)` returns the variance of the statistics `stats`.
__attribute__ ((__nonnull__, __pure__))
double stats_var(struct stats const*);

// The call `stats_sd(stats)` returns
// the standard deviation of the statistics `stats`.
__attribute__ ((__nonnull__, __pure__))
double stats_sd(struct stats const*);

// The call `stats_sem(stats)` returns
// the standard error of the mean of the statistics `stats`.
__attribute__ ((__nonnull__, __pure__))
double stats_sem(struct stats const*);

#endif
