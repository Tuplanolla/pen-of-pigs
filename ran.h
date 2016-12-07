#ifndef RAN_H
#define RAN_H

#include "exts.h"
#include <gsl/gsl_rng.h>
#include <stddef.h>

// The call `ran_index(rng, n)` uses `rng` to generate
// a random index for an array of size `n`.
__attribute__ ((__nonnull__))
size_t ran_index(gsl_rng*, size_t);

// The call `ran_uopen(rng, x)` uses `rng` to generate
// a random floating-point number in the open interval from `0` to `x`.
__attribute__ ((__nonnull__))
double ran_uopen(gsl_rng*, double);

// The call `ran_open(rng, x)` uses `rng` to generate
// a random floating-point number in the open interval from `-x` to `x`.
__attribute__ ((__nonnull__))
double ran_open(gsl_rng*, double);

// The call `ran_dateid(rng, pre, buf, n)` uses `rng` to generate
// an identifier consisting of the prefix `pre`,
// the current date and a random hexadecimal value,
// writes it into the buffer `buf` of size `n` and returns `true`.
// Otherwise `false` is returned.
__attribute__ ((__nonnull__))
bool ran_dateid(gsl_rng*, char const*, char*, size_t);

#endif
