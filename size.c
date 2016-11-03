#include "size.h"
#include <stddef.h>

int size_cmp(size_t const x, size_t const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

size_t size_min(size_t const x, size_t const y) {
  return x < y ? x : y;
}

size_t size_max(size_t const x, size_t const y) {
  return x > y ? x : y;
}

size_t size_pow(size_t const x, size_t const y) {
  size_t z = 1;

  for (size_t i = 0;
      i < y;
      ++i)
    z *= x;

  return z;
}
