#include "exts.h"
#include "fp.h"
#include <math.h>

double fp_identity(double const x) {
  return x;
}

double fp_constant(double const x,
    __attribute__ ((__unused__)) double const y) {
  return x;
}

double fp_zero(__attribute__ ((__unused__)) double const x) {
  return 0.0;
}

double fp_one(__attribute__ ((__unused__)) double const x) {
  return 1.0;
}

int fp_cmp(double const x, double const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

double fp_midpoint(double const x, double const y) {
  return (x + y) / 2.0;
}

fp_div_t fp_div(double const x, double const y) {
  fp_div_t const z = {
    .quot = trunc(x / y),
    .rem = fmod(x, y)
  };

  return z;
}

double fp_min(double const x, double const y) {
  return fmin(x, y);
}

double fp_max(double const x, double const y) {
  return fmax(x, y);
}

double fp_pow(double const x, double const y) {
  return pow(x, y);
}

double fp_uwrap(double const x, double const y) {
  return x - y * floor(x / y);
}

double fp_wrap(double const x, double const y) {
  return x - y * nearbyint(x / y);
}

double fp_lerp(double const x0, double const x1,
    double const y0, double const y1,
    double const x) {
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}
