#include "exts.h"
#include "flt.h"
#include <math.h>

int cmp(double const x, double const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

double wrap(double const x, double const y) {
  return fmod(y + fmod(x, y), y);
}

double midpoint(double const x, double const y) {
  return (x + y) / 2;
}

double lerp(double const x0, double const x1,
    double const y0, double const y1,
    double const x) {
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

double identity(double const x) {
  return x;
}

double zero(__attribute__ ((__unused__)) double const x) {
  return 0;
}

double one(__attribute__ ((__unused__)) double const x) {
  return 1;
}

double periodic(double const x, double const p) {
  double const p2 = p / 2;
  return x > p2 ? x - p : x < -p2 ? x + p : x;
}
