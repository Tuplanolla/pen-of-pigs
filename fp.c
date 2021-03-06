#include "exts.h"
#include "fp.h"
#include <gsl/gsl_math.h>
#include <math.h>
#include <stddef.h>

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
  fp_div_t const qr = {
    .quot = trunc(x / y),
    .rem = fmod(x, y)
  };

  return qr;
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

double fp_log(double const x, double const y) {
  return log2(x) / log2(y);
}

double fp_rt(double const x, double const y) {
  return pow(x, 1.0 / y);
}

double fp_clamp(double const x, double const a, double const b) {
  return x < a ? a : x >= b ? b : x;
}

double fp_sclamp(double const x, double const b) {
  double const a = b / 2.0;

  return x < -a ? a : x >= a ? a : x;
}

double fp_uclamp(double const x, double const b) {
  return x < 0.0 ? 0.0: x >= b ? b : x;
}

double fp_wrap(double const x, double const a, double const b) {
  double const c = b - a;

  return x - c * floor((x - a) / c);
}

double fp_swrap(double const x, double const b) {
  return x - b * nearbyint(x / b);
}

double fp_uwrap(double const x, double const b) {
  return x - b * floor(x / b);
}

double fp_lerp(double const x,
    double const x0, double const x1, double const y0, double const y1) {
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

double fp_lorp(double const x,
    double const x0, double const x1, double const y0, double const y1) {
  return log(fp_lerp(exp(x), exp(x0), exp(x1), exp(y0), exp(y1)));
}

// The recurrence relation
// $V_0(r) = 1$, $V_1(r) = 2 r$, $V_n(r) = (\twopi r^2 / n) V_{n - 2}(r)$
// is applied.
double fp_ball_volume(double const r, size_t d) {
  // return d == 0 ? 1.0 : d == 1 ? 2.0 * r :
  //   (M_2PI * gsl_pow_2(r) / (double) d) * fp_ball_volume(r, d - 2);

  double v = d % 2 == 0 ? 1.0 : 2.0 * r;

  while (d > 1) {
    v *= M_2PI * gsl_pow_2(r) / (double) d;
    d -= 2;
  }

  return v;
}

// The recurrence relation
// $A_0(r) = 0$, $A_1(r) = 2$, $A_2(r) = \twopi r$,
// $A_n(r) = (\twopi r^2 / (n - 2)) A_{n - 2}(r)$
// is applied.
double fp_ball_surfarea(double const r, size_t d) {
  // return d == 0 ? 0.0 : d == 1 ? 2.0 : d == 2 ? M_2PI * r :
  //   (M_2PI * gsl_pow_2(r) / (double) (d - 2)) * fp_ball_surfarea(r, d - 2);

  double a = d == 0 ? 0.0 : d % 2 == 1 ? 2.0 : M_2PI * r;

  while (d > 2) {
    d -= 2;
    a *= M_2PI * gsl_pow_2(r) / (double) d;
  }

  return a;
}

double fp_dbalance(double const a, double const r) {
  return 0.5 + a / (a + r);
}

double fp_cbalance(double const a, double const r) {
  return 1.0 - exp(-a) + exp(-r);
}
