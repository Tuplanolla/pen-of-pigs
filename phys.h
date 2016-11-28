#ifndef PHYS_H
#define PHYS_H

#include "exts.h"

extern double const hbar;
extern double const kB;

// Lennard--Jones potential 6--12.
__attribute__ ((__const__))
double lj612(double, double, double);

// Lennard--Jones 6--12 potential from squared arguments (for speed).
__attribute__ ((__const__))
double lj6122(double, double, double);

// Harmonic potential 6--12.
__attribute__ ((__const__))
double harm(double, double, double);

// Harmonic potential 6--12 from squared arguments (for speed).
__attribute__ ((__const__))
double harm2(double, double, double);

// The call `hermite(n, x)` returns the value
// of the `n`th Hermite polynomial at `x`.
// Note that the function is only defined when `0 <= n < 8`,
// because the programmer is lazy.
__attribute__ ((__const__))
double hermite(unsigned int, double);

#endif
