#ifndef SIZE_H
#define SIZE_H

#include "exts.h"
#include <limits.h>
#include <stddef.h>

#define SIZE_MAX_BIT (CHAR_BIT * sizeof (size_t))
#define SQRT_SIZE_MAX_BIT (SIZE_MAX_BIT / 2)
#define SQRT_SIZE_MAX ((size_t) 1 << SQRT_SIZE_MAX_BIT)

// The statement `n = size_identity(n)` does not do anything.
// This is analogous to `fp_identity`.
__attribute__ ((__const__, __pure__))
size_t size_identity(size_t);

// The call `size_zero(n)` returns `0`.
// This is analogous to `fp_zero`.
__attribute__ ((__const__, __pure__))
size_t size_zero(size_t);

// The call `size_one(n)` returns `1`.
// This is analogous to `fp_one`.
__attribute__ ((__const__, __pure__))
size_t size_one(size_t);

// The call `size_cmp(n, k)` returns
//
// * `-1` if `n < k`,
// * `1` if `n > k` and
// * `0` otherwise.
//
// This is analogous to `fp_cmp`.
__attribute__ ((__const__, __pure__))
int size_cmp(size_t, size_t);

// The call `size_midpoint(n, k)` returns the arithmetic mean of `n` and `k`.
// This is analogous to `fp_midpoint`.
__attribute__ ((__const__, __pure__))
size_t size_midpoint(size_t, size_t);

// This structure holds the quotient and remainder of a division
// in unspecified order.
typedef struct {
  size_t quot;
  size_t rem;
} size_div_t;

// The statement `m = size_div(n, k)` solves
// the division equation `m.quot * k + m.rem == n` for `m`,
// where `m.quot` is the quotient and `m.rem` is the remainder
// of the expression `n / k`.
// This is analogous to `div`.
__attribute__ ((__const__, __pure__))
size_div_t size_div(size_t, size_t);

// The call `size_min(n, k)` returns the lesser of `n` and `k`.
// This is analogous to `fmin`.
__attribute__ ((__const__, __pure__))
size_t size_min(size_t, size_t);

// The call `size_max(n, k)` returns the greater of `n` and `k`.
// This is analogous to `fmax`.
__attribute__ ((__const__, __pure__))
size_t size_max(size_t, size_t);

// The call `size_pow(n, k)` returns `n` raised to the power of `k`.
// This is analogous to `pow`.
__attribute__ ((__const__, __pure__))
size_t size_pow(size_t, size_t);

// The call `size_filog(n, k)` returns the floor
// of the base `k` logarithm of `n`.
// This is analogous to `fp_log`.
__attribute__ ((__const__, __pure__))
size_t size_filog(size_t, size_t);

// The call `size_cilog(n, k)` returns the ceiling
// of the base `k` logarithm of `n`.
// This is analogous to `fp_log`.
__attribute__ ((__const__, __pure__))
size_t size_cilog(size_t, size_t);

// The call `size_firt(n, k)` returns the floor of the `k`th root of `n`.
// This is analogous to `fp_rt`.
// Note that the result may be wrong for large arguments.
__attribute__ ((__const__, __pure__))
size_t size_firt(size_t, size_t);

// The call `size_cirt(n, k)` returns the ceiling of the `k`th root of `n`.
// This is analogous to `fp_rt`.
// Note that the result may be wrong for large arguments.
__attribute__ ((__const__, __pure__))
size_t size_cirt(size_t, size_t);

// The call `size_hc(ilin, nlin, icub, ndim)` sets `icub`
// to the hypercubical index corresponding to the linear index `ilin`
// in a hypercube of side length `nlin` and dimension `ndim`.
// Immediately after the call it is guaranteed that
// `ilin == size_unhc(nlin, icub, ndim)`.
__attribute__ ((__nonnull__))
void size_hc(size_t, size_t, size_t*, size_t);

// The statement `ilin = size_unhc(nlin, icub, ndim)` sets `ilin`
// to the linear index corresponding to the hypercubical index `icub`
// in a hypercube of side length `nlin` and dimension `ndim`.
// Immediately after the statement it is guaranteed that
// `size_hc(ilin, nlin, icub, ndim)` does nothing.
__attribute__ ((__nonnull__, __pure__))
size_t size_unhc(size_t, size_t const*, size_t);

// The call `size_gcd(n, k)` returns the greatest common divisor of
// `n` and `k` or `0` if both `n` and `k` are zero.
__attribute__ ((__const__, __pure__))
size_t size_gcd(size_t, size_t);

// The call `size_uclamp(n, k)` returns
//
// * `n` if `0 <= n < k` and
// * `k - 1` if `n >= k`.
//
// This is analogous to `fp_uclamp`.
__attribute__ ((__const__, __pure__))
size_t size_uclamp(size_t, size_t);

// The statement `m = size_uwrap(n, k)` solves
// the periodic equation `m == n - p * k` for `m`,
// where `0 <= m < k` and `p` is some integer.
// This is analogous to `fp_uwrap`.
__attribute__ ((__const__, __pure__))
size_t size_uwrap(size_t, size_t);

// The statement `m = size_uwrap_inc(n, k)` is equivalent
// to `m = size_uwrap(n + 1, k)` without overflowing (or wrapping).
__attribute__ ((__const__, __pure__))
size_t size_uwrap_inc(size_t, size_t);

// The statement `m = size_uwrap_dec(n, k)` is equivalent
// to `m = size_uwrap(n - 1, k)` without underflowing (or wrapping).
__attribute__ ((__const__, __pure__))
size_t size_uwrap_dec(size_t, size_t);

#endif
