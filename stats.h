#ifndef STATS_H
#define STATS_H

#include "exts.h"
#include <stddef.h>

struct stats;

void stats_free(struct stats*);

__attribute__ ((__malloc__))
struct stats* stats_alloc(void);

__attribute__ ((__nonnull__))
void stats_accum(struct stats*, double);

__attribute__ ((__nonnull__, __pure__))
double stats_mean(struct stats const*);

__attribute__ ((__nonnull__, __pure__))
double stats_var(struct stats const*);

__attribute__ ((__nonnull__, __pure__))
double stats_sd(struct stats const*);

__attribute__ ((__nonnull__, __pure__))
double stats_sem(struct stats const*);

#endif
