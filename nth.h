#ifndef NTH_H
#define NTH_H

#include "exts.h"
#include <stddef.h>

/*
The call `nth_prime(n)` returns
the `n`th prime number if `1 <= n < m` where
`m` is the smallest guaranteed `UINT_MAX` or
`1` if `n == 0`.
*/
__attribute__ ((__const__))
unsigned int nth_prime(size_t n);

#endif
