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
The call `z = size_wrap(x, y)` returns such `z` that
`0 <= z < y` and `z = x - n * y` for some integer `n`.
*/
__attribute__ ((__const__))
size_t size_wrap(size_t x, size_t y);

/*
The call `size_midpoint(x, y)` returns `(x + y) / 2`.
*/
__attribute__ ((__const__))
size_t size_midpoint(size_t x, size_t y);

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

typedef struct {
  size_t quot;
  size_t rem;
} size_div_t;

/*
The call `z = size_div(x, y)` returns
the quotient `z.quot` and the remainder `z.rem`
of the expression `x / y` in unspecified order.
*/
__attribute__ ((__const__))
size_div_t size_div(size_t x, size_t y);

/*
The call `size_pow(x, y)` returns `x` to the power of `y`.
*/
__attribute__ ((__const__))
size_t size_pow(size_t x, size_t y);

#endif
