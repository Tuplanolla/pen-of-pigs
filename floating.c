#include "floating.h"
#include "report.h"
#include <gsl/gsl_math.h>

int cmp(double const x, double const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

double lerp(double const x0, double const x1,
    double const y0, double const y1,
    double const x) {
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

double hermite(unsigned int const n, double const x) {
  switch (n) {
    case 0:
      return 1;
    case 1:
      return 2 * x;
    case 2:
      return 4 * gsl_pow_uint(x, 2)
        - 2;
    case 3:
      return 8 * gsl_pow_uint(x, 3)
        - 12 * x;
    case 4:
      return 16 * gsl_pow_uint(x, 4)
        - 48 * gsl_pow_uint(x, 2)
        + 12;
    case 5:
      return 32 * gsl_pow_uint(x, 5)
        - 160 * gsl_pow_uint(x, 3)
        + 120 * x;
    case 6:
      return 64 * gsl_pow_uint(x, 6)
        - 480 * gsl_pow_uint(x, 4)
        + 720 * gsl_pow_uint(x, 2)
        - 120;
    case 7:
      return 128 * gsl_pow_uint(x, 7)
        - 1344 * gsl_pow_uint(x, 5)
        + 3360 * gsl_pow_uint(x, 3)
        - 1680 * x;
    default:
      /* return 2 * x * hermite(n - 1, x) - dhermite(n - 1, x); */
      warn(hermite);
      return NAN;
  }
}
