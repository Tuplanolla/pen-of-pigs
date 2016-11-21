#ifndef FP_H
#define FP_H

#include "exts.h"

#define M_2PI 6.283185307179586

/*
The call `fp_cmp(x, y)` returns

* `-1` if `x < y`,
* `1` if `x > y` and
* `0` otherwise.
*/
__attribute__ ((__const__))
int fp_cmp(double, double);

/*
The statement `z = fp_wrap(x, y)` solves
the fp_periodic equation `z == x - n * y` for `z`,
where `0 <= z < y` and `n` is some integer.
*/
__attribute__ ((__const__))
double fp_wrap(double, double);

/*
The call `fp_midpoint(x, y)` returns the arithmetic mean of `x` and `y`.
*/
__attribute__ ((__const__))
double fp_midpoint(double, double);

/*
The statement `y = fp_lerp(x0, x1, y0, y1, x)` solves
the linear interpolation equation
`(x - x0) / (x1 - x0) == (y - y0) / (y1 - y0)` for `y`.
*/
__attribute__ ((__const__))
double fp_lerp(double, double, double, double, double);

/*
The statement `x = fp_identity(x)` does not do anything.
*/
__attribute__ ((__const__))
double fp_identity(double);

/*
The call `fp_zero(x)` returns `0`.
*/
__attribute__ ((__const__))
double fp_zero(double);

/*
The call `fp_one(x)` returns `1`.
*/
__attribute__ ((__const__))
double fp_one(double);

/*
The call `fp_periodic(x, p)` minimizes the periodicity expression `x + n * p`,
where `n` is some integer.
*/
__attribute__ ((__const__))
double fp_periodic(double, double);

#endif
