#ifndef FLT_H
#define FLT_H

#include "exts.h"

#define M_2PI 6.283185307179586

/*
The call `cmp(x, y)` returns

* `-1` if `x < y`,
* `1` if `x > y` and
* `0` otherwise.
*/
__attribute__ ((__const__))
int cmp(double, double);

/*
The statement `z = wrap(x, y)` solves
the periodic equation `z == x - n * y` for `z`,
where `0 <= z < y` and `n` is some integer.
*/
__attribute__ ((__const__))
double wrap(double, double);

/*
The call `midpoint(x, y)` returns the arithmetic mean of `x` and `y`.
*/
__attribute__ ((__const__))
double midpoint(double, double);

/*
The statement `y = lerp(x0, x1, y0, y1, x)` solves
the linear interpolation equation
`(x - x0) / (x1 - x0) == (y - y0) / (y1 - y0)` for `y`.
*/
__attribute__ ((__const__))
double lerp(double, double, double, double, double);

/*
The statement `x = identity(x)` does not do anything.
*/
__attribute__ ((__const__))
double identity(double);

/*
The call `zero(x)` returns `0`.
*/
__attribute__ ((__const__))
double zero(double);

/*
The call `one(x)` returns `1`.
*/
__attribute__ ((__const__))
double one(double);

/*
The statements `y = periodic(x, p)` sets `y` to `x` with period `p`.
This description is too vague.
*/
__attribute__ ((__const__))
double periodic(double, double);

#endif
