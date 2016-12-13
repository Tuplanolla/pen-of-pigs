#include "opt.h"
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef _GNU_SOURCE

#include <getopt.h>

int opt_parse(int const argc, char* const* const argv,
    char const* const shortstr, char const* const* const longstr) {
  if (longstr == NULL)
    return getopt(argc, argv, shortstr);
  else {
    size_t m = 0;

    for (size_t i = 0; shortstr[i] != '\0'; ++i)
      if (shortstr[i] != ':')
        ++m;

    struct option* const longopt = malloc(m * sizeof *longopt);

    for (size_t i = 0, j = 0; j < m; ++j) {
      size_t k = 0;

      while (shortstr[i + k + 1] == ':')
        ++k;

      longopt[j].name = longstr[j];
      longopt[j].has_arg = (int) k;
      longopt[j].flag = NULL;
      longopt[j].val = shortstr[i];

      i += k + 1;
    }

    int const i = getopt_long(argc, argv, shortstr, longopt, NULL);

    free(longopt);

    return i;
  }
}

#else

int opt_parse(int const argc, char* const* const argv,
    char const* const shortstr,
    __attribute__ ((__unused__)) char const* const* const longstr) {
  return getopt(argc, argv, shortstr);
}

#endif

bool opt_parse_size(size_t* const n, size_t const min, size_t const max) {
  char* endptr;
  errno = 0;
  unsigned long long int const k = strtoull(optarg, &endptr, 10);
  if (errno != 0 || endptr == optarg || *endptr != '\0') {
    errno = EINVAL;

    return false;
  } else if (k > SIZE_MAX || (size_t) k < min || (size_t) k > max) {
    errno = ERANGE;

    return false;
  } else
    *n = (size_t) k;

  return true;
}

bool opt_parse_fp(double* const x, double const min, double const max) {
  char* endptr;
  errno = 0;
  double const y = strtod(optarg, &endptr);
  if (errno != 0 || endptr == optarg || *endptr != '\0') {
    errno = EINVAL;

    return false;
  } else if (y < min || y > max) {
    errno = ERANGE;

    return false;
  } else
    *x = y;

  return true;
}
