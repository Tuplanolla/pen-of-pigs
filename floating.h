#ifndef FLOATING_H
#define FLOATING_H

#include "exts.h"

#define M_2PI 6.283185307179586

/*
The call `size_cmp(x, y)` returns

* `-1` if `x < y`,
* `1` if `x > y` and
* `0` otherwise.
*/
__attribute__ ((__const__))
int cmp(double x, double y);

/*
The call `lerp(x0, x1, y0, y1, x)` returns the solution
for `y` in the linear interpolation equation
`(x - x0) / (x1 - x0) == (y - y0) / (y1 - y0)`.
*/
__attribute__ ((__const__))
double lerp(double x0, double x1, double y0, double y1, double x);

/*
The call `hermite(n, x)` returns the value
of the `n`th Hermite polynomial at `x`.
Note that the function is only defined when `0 <= n < 8`.
*/
__attribute__ ((__const__))
double hermite(unsigned int n, double x);

#endif
