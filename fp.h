#ifndef FP_H
#define FP_H

#include "exts.h"

#define M_2PI 6.283185307179586

// The statement `x = fp_identity(x)` does not do anything.
// This is analogous to `size_identity`.
__attribute__ ((__const__))
double fp_identity(double);

// The statement `x = fp_constant(x, y)` does not do anything.
// This is analogous to `size_constant`.
__attribute__ ((__const__))
double fp_constant(double, double);

// The call `fp_zero(x)` returns `0`.
// This is analogous to `size_zero`.
__attribute__ ((__const__))
double fp_zero(double);

// The call `fp_one(x)` returns `1`.
// This is analogous to `size_one`.
__attribute__ ((__const__))
double fp_one(double);

// The call `fp_cmp(x, y)` returns
//
// * `-1` if `x < y`,
// * `1` if `x > y` and
// * `0` otherwise.
//
// This is analogous to `size_cmp`.
__attribute__ ((__const__))
int fp_cmp(double, double);

// The call `fp_midpoint(x, y)` returns the arithmetic mean of `x` and `y`.
// This is analogous to `size_midpoint`.
__attribute__ ((__const__))
double fp_midpoint(double, double);

// This structure holds the quotient and remainder of a division
// in unspecified order.
typedef struct {
  double quot;
  double rem;
} fp_div_t;

// The statement `z = fp_div(x, y)` solves
// the division equation `z.quot * y + z.rem == x` for `z`,
// where `z.quot` is the quotient and `z.rem` is the remainder
// of the expression `x / y`.
// This is analogous to `div` or `size_div`.
__attribute__ ((__const__))
fp_div_t fp_div(double, double);

// The statement `fp_min(x, y)` returns the lesser of `x` and `y`.
// This is analogous to `fmin` or `size_min`.
__attribute__ ((__const__))
double fp_min(double, double);

// The statement `fp_max(x, y)` returns the greater of `x` and `y`.
// This is analogous to `fmax` or `size_max`.
__attribute__ ((__const__))
double fp_max(double, double);

// The call `fp_pow(x, y)` returns `x` raised to the power of `y`.
// This is analogous to `pow` or `size_pow`.
__attribute__ ((__const__))
double fp_pow(double, double);

// The call `fp_log(x, y)` returns the the base `y` logarithm of `x`.
// This is analogous to `size_filog` or `size_cilog`.
__attribute__ ((__const__))
double fp_log(double, double);

// The call `fp_rt(x, y)` returns the `y`th root of `x`.
// This is analogous to `size_firt` or `size_cirt`.
__attribute__ ((__const__))
double fp_rt(double, double);

// The call `fp_clamp(x, a, b)` returns
//
// * `x` if `a <= x < b`,
// * `a` if `x < a` and
// * `b` if `x >= b`.
//
// This is equivalent to `fp_uclamp(x, b)` when `a == 0`.
__attribute__ ((__const__))
double fp_clamp(double, double, double);

// The call `fp_uclamp(x, b)` returns
//
// * `x` if `0 <= x < b`,
// * `0` if `x < 0` and
// * `b` if `x > b`.
//
// This is analogous to `size_uclamp`.
__attribute__ ((__const__))
double fp_uclamp(double, double);

// The statement `z = fp_wrap(x, y)` solves
// the periodic equation `z == x + n * y` for `z`,
// where `-y / 2 <= z < y / 2` and `n` is some integer.
__attribute__ ((__const__))
double fp_wrap(double, double);

// The statement `z = fp_uwrap(x, y)` solves
// the periodic equation `z == x + n * y` for `z`,
// where `0 <= z < y` and `n` is some integer.
// This is analogous to `size_uwrap`.
__attribute__ ((__const__))
double fp_uwrap(double, double);

// The statement `y = fp_lerp(x, x0, x1, y0, y1)` solves
// the linear interpolation equation
// `(x - x0) / (x1 - x0) == (y - y0) / (y1 - y0)` for `y`.
__attribute__ ((__const__))
double fp_lerp(double, double, double, double, double);

// The statement `y = fp_lorp(x, x0, x1, y0, y1)` is equivalent to
// `y = log(fp_lerp(exp(x), exp(x0), exp(x1), exp(y0), exp(y1)))`.
__attribute__ ((__const__))
double fp_lorp(double, double, double, double, double);

// The statement `b = fp_dbalance(a, r)` sets `b` such that
//
// * `b > 0 && b < 1` if `a > r`,
// * `b > 1` if `a < r` and
// * `b == 1` otherwise.
//
// The function `fp_dbalance` is essentially
// smooth and continuous everywhere except for the origin.
// Iterating `fp_dbalance` is faster than `fp_cbalance`,
// but tends to diverge into bounded oscillations.
__attribute__ ((__const__))
double fp_dbalance(double, double);

// The statement `b = fp_cbalance(a, r)` sets `b` such that
//
// * `b > 0 && b < 1` if `a > r`,
// * `b > 1` if `a < r` and
// * `b == 1` otherwise.
//
// The function `fp_cbalance` is essentially
// smooth and continuous everywhere.
// Iterating `fp_cbalance` is slower than `fp_dbalance`,
// but tends to converge quickly.
__attribute__ ((__const__))
double fp_cbalance(double, double);

#endif
