#ifndef EST_H
#define EST_H

#include "exts.h"
#include <stddef.h>

struct est;

void est_free(struct est*);

__attribute__ ((__malloc__))
struct est* est_alloc(void);

__attribute__ ((__nonnull__))
void est_accum(struct est*, double);

__attribute__ ((__nonnull__, __pure__))
double est_mean(struct est const*);

__attribute__ ((__nonnull__, __pure__))
double est_var(struct est const*);

__attribute__ ((__nonnull__, __pure__))
double est_sd(struct est const*);

__attribute__ ((__nonnull__, __pure__))
double est_sem(struct est const*);

#endif
