#include "exts.h"
#include "size.h"
#include <stddef.h>

size_t size_identity(size_t const n) {
  return n;
}

size_t size_zero(__attribute__ ((__unused__)) size_t const n) {
  return 0;
}

size_t size_one(__attribute__ ((__unused__)) size_t const n) {
  return 1;
}

int size_cmp(size_t const n, size_t const k) {
  return n < k ? -1 : n > k ? 1 : 0;
}

size_t size_midpoint(size_t const n, size_t const k) {
  // Note that `(n + k) / 2` could overflow (or wrap) and
  // `n < k ? n + (k - n) / 2 : k + (n - k) / 2` could have bad performance.
  return n / 2 + k / 2 + (n % 2 + k % 2) / 2;
}

size_div_t size_div(size_t const n, size_t const k) {
  size_div_t const qr = {
    .quot = n / k,
    .rem = n % k
  };

  return qr;
}

size_t size_min(size_t const n, size_t const k) {
  return n < k ? n : k;
}

size_t size_max(size_t const n, size_t const k) {
  return n > k ? n : k;
}

size_t size_pow(size_t const n, size_t const k) {
  size_t m = 1;

  for (size_t i = 0; i < k; ++i)
    m *= n;

  return m;
}

size_t size_filog(size_t n, size_t const k) {
  dynamic_assert(n <= 0, "invalid argument");
  dynamic_assert(k <= 1, "invalid base");

  size_t m = 0;

  while (n >= k) {
    n /= k;
    ++m;
  }

  return m;
}

size_t size_cilog(size_t const n, size_t const k) {
  dynamic_assert(n <= 0, "invalid argument");
  dynamic_assert(k <= 1, "invalid base");

  return n <= 1 ? 0 : size_filog(n - 1, k) + 1;
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

void size_hc(size_t ilin, size_t const nlin,
    size_t* const icub, size_t const ndim) {
  for (size_t idim = 0; idim < ndim; ++idim) {
    size_div_t const qr = size_div(ilin, nlin);

    ilin = qr.quot;
    icub[idim] = qr.rem;
  }
}

size_t size_unhc(size_t const nlin,
    size_t const* const icub, size_t const ndim) {
  size_t ilin = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    ilin *= nlin;
    ilin += icub[ndim - 1 - idim];
  }

  return ilin;
}

size_t size_uclamp(size_t const n, size_t const k) {
  return n >= k ? k - 1 : n;
}

size_t size_uwrap(size_t const n, size_t const k) {
  return n % k;
}

size_t size_uwrap_inc(size_t const n, size_t const k) {
  return n == k - 1 ? 0 : n + 1;
}

size_t size_uwrap_dec(size_t const n, size_t const k) {
  return n == 0 ? k - 1 : n - 1;
}
