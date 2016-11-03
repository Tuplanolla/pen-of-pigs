#ifndef SIZE_H
#define SIZE_H

#include "exts.h"
#include <stddef.h>

/*
The call `size_cmp(x, y)` returns

* `-1` if `x < y`,
* `1` if `x > y` and
* `0` otherwise.
*/
__attribute__ ((__const__))
int size_cmp(size_t x, size_t y);

/*
The call `size_min(x, y)` returns the lesser of `x` and `y`.
*/
__attribute__ ((__const__))
size_t size_min(size_t x, size_t y);

/*
The call `size_max(x, y)` returns the greater of `x` and `y`.
*/
__attribute__ ((__const__))
size_t size_max(size_t x, size_t y);

/*
The call `size_pow(x, y)` returns `x` to the power of `y`.
*/
__attribute__ ((__const__))
size_t size_pow(size_t x, size_t y);

#endif
