#ifndef SIZE_H
#define SIZE_H

#include "exts.h"
#include <limits.h>
#include <stddef.h>

#define RT_SIZE_MAX(n) ((size_t) 1 << CHAR_BIT * sizeof (size_t) / (n))

// The statement `x = size_identity(x)` does not do anything.
// This is analogous to `fp_identity`.
__attribute__ ((__const__))
size_t size_identity(size_t);

// The call `size_zero(x)` returns `0`.
// This is analogous to `fp_zero`.
__attribute__ ((__const__))
size_t size_zero(size_t);

// The call `y = size_one(x)` returns `1`.
// This is analogous to `fp_one`.
__attribute__ ((__const__))
size_t size_one(size_t);

// The call `size_cmp(x, y)` returns
//
// * `-1` if `x < y`,
// * `1` if `x > y` and
// * `0` otherwise.
//
// This is analogous to `fp_cmp`.
__attribute__ ((__const__))
int size_cmp(size_t, size_t);

// The call `size_midpoint(x, y)` returns the arithmetic mean of `x` and `y`.
// This is analogous to `fp_midpoint`.
__attribute__ ((__const__))
size_t size_midpoint(size_t, size_t);

// This structure holds the quotient and remainder of a division
// in unspecified order.
typedef struct {
  size_t quot;
  size_t rem;
} size_div_t;

// The statement `z = size_div(x, y)` solves
// the division equation `z.quot * y + z.rem == x` for `z`,
// where `z.quot` is the quotient and `z.rem` is the remainder
// of the expression `x / y`.
// This is analogous to `div`.
__attribute__ ((__const__))
size_div_t size_div(size_t, size_t);

// The statement `size_min(x, y)` returns the lesser of `x` and `y`.
// This is analogous to `fmin`.
__attribute__ ((__const__))
size_t size_min(size_t, size_t);

// The statement `size_max(x, y)` returns the greater of `x` and `y`.
// This is analogous to `fmax`.
__attribute__ ((__const__))
size_t size_max(size_t, size_t);

// The call `size_pow(x, y)` returns `x` raised to the power of `y`.
// This is analogous to `pow`.
__attribute__ ((__const__))
size_t size_pow(size_t, size_t);

// The statement `z = size_uwrap(x, y)` solves
// the periodic equation `z == x - n * y` for `z`,
// where `0 <= z < y` and `n` is some integer.
// This is analogous to `fp_uwrap`.
__attribute__ ((__const__))
size_t size_uwrap(size_t, size_t);

#endif
