#ifndef ALLOC_H
#define ALLOC_H

#include "exts.h"
#include <stddef.h>

typedef void* (* alloc_proc)(size_t);

// The statement `p = alloc_err(n)` sets `p` to `NULL`.
// This is analogous to `malloc` that always fails.
__attribute__ ((__malloc__))
void* alloc_err(size_t);

// The call `alloc_ran_rel(q)` sets the reliability of `alloc_ran` to `q`.
void alloc_ran_rel(double);

// The statement `p = alloc_ran(n)` sets the pointer `p`
// to a newly allocated memory block of `n` bytes or
// fails artificially and sets `p` to `NULL`.
// The reliability `q` or the equivalent inverse probability of failure `1 - q`
// can be set with the call `alloc_ran_rel(q)`.
// This is analogous to `malloc` that sometimes fails.
__attribute__ ((__malloc__))
void* alloc_ran(size_t);

#endif
