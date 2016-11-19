#ifndef NTH_H
#define NTH_H

#include "exts.h"
#include <stddef.h>

/*
The call `p = nth_prime(n)` returns the `(n + 1)`th prime number
if `0 <= p < m` and `m` is the smallest guaranteed `UINT_MAX`.
*/
__attribute__ ((__const__))
unsigned int nth_prime(size_t);

#endif
