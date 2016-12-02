#ifndef RAN_H
#define RAN_H

#include "exts.h"
#include <gsl/gsl_rng.h>
#include <stddef.h>

size_t ran_index(gsl_rng*, size_t);

__attribute__ ((__nonnull__))
double ran_uopen(gsl_rng*, double);

__attribute__ ((__nonnull__))
double ran_open(gsl_rng*, double);

#endif
