#ifndef FLT_H
#define FLT_H

#include "exts.h"

#define M_2PI 6.283185307179586

/*
The call `size_cmp(x, y)` returns

* `-1` if `x < y`,
* `1` if `x > y` and
* `0` otherwise.
*/
__attribute__ ((__const__))
int cmp(double, double);

/*
The call `wrap(x, y)` returns such `z` that
`0 <= z < y` and `z = x - n * y` for some integer `n`.
*/
__attribute__ ((__const__))
double wrap(double, double);

/*
The call `midpoint(x, y)` returns the mean of `x` and `y`.
*/
__attribute__ ((__const__))
double midpoint(double, double);

/*
The call `lerp(x0, x1, y0, y1, x)` returns the solution
for `y` in the linear interpolation equation
`(x - x0) / (x1 - x0) == (y - y0) / (y1 - y0)`.
*/
__attribute__ ((__const__))
double lerp(double, double, double, double, double);

/*
The call `hermite(n, x)` returns the value
of the `n`th Hermite polynomial at `x`.
Note that the function is only defined when `0 <= n < 8`.
*/
__attribute__ ((__const__))
double hermite(unsigned int, double);

#endif
