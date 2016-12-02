#include "exts.h"
#include "size.h"
#include <stddef.h>

size_t size_identity(size_t const x) {
  return x;
}

size_t size_zero(__attribute__ ((__unused__)) size_t const x) {
  return 0;
}

size_t size_one(__attribute__ ((__unused__)) size_t const x) {
  return 1;
}

int size_cmp(size_t const x, size_t const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

size_t size_midpoint(size_t const x, size_t const y) {
  // Note that `(x + y) / 2` could overflow (or wrap) and
  // `x < y ? x + (y - x) / 2 : y + (x - y) / 2` could have bad performance.
  return x / 2 + y / 2 + (x % 2 + y % 2) / 2;
}

size_div_t size_div(size_t const x, size_t const y) {
  size_div_t const z = {
    .quot = x / y,
    .rem = x % y
  };

  return z;
}

size_t size_min(size_t const x, size_t const y) {
  return x < y ? x : y;
}

size_t size_max(size_t const x, size_t const y) {
  return x > y ? x : y;
}

size_t size_pow(size_t const x, size_t const y) {
  size_t z = 1;

  for (size_t i = 0; i < y; ++i)
    z *= x;

  return z;
}

size_t size_firt(size_t const n, size_t const k) {
  if (n <= 1)
    return n;
  else {
    size_t const p = k - 1;
    size_t r = n + 1;
    size_t m = n;

    while (m < r) {
      r = m;
      m = (p * r + n / size_pow(r, p)) / k;
    }

    return r;
  }
}

size_t size_cirt(size_t const n, size_t const k) {
  return n <= 1 ? n : size_firt(n - 1, k) + 1;
}

size_t size_uwrap(size_t const x, size_t const y) {
  return x % y;
}

size_t size_uwrap_inc(size_t const x, size_t const y) {
  return x == y - 1 ? 0 : x + 1;
}

size_t size_uwrap_dec(size_t const x, size_t const y) {
  return x == 0 ? y - 1 : x - 1;
}
