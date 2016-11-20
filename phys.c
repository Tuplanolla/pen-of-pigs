#include "err.h"
#include "phys.h"
#include <gsl/gsl_math.h>

double lj612(double const r, double const sigma, double const epsilon) {
  double const sigmar6 = gsl_pow_6(sigma / r);

  return 4 * epsilon * (gsl_pow_2(sigmar6) - sigmar6);
}

double lj6122(double const r2, double const sigma2, double const epsilon) {
  double const sigmar2 = sigma2 / r2;

  return 4 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

double harm(double const r, double const omega, double const m) {
  return m * gsl_pow_2(omega * r) / 2;
}

double harm2(double const r2, double const omega2, double const m) {
  return m * omega2 * r2 / 2;
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
      err_abort(hermite);
  }
}
