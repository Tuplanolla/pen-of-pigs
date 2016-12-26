#ifndef OPT_H
#define OPT_H

#include "exts.h"
#include <stdbool.h>
#include <stddef.h>

// The statement `i = opt_parse(argc, argv, shortstr, longstr)` parses
// one argument from the command-line argument array `argv` of length `argc`,
// where `shortstr` contains the short options
// in a format understood by `getopt` and
// `longstr` either contains the corresponding long options or is `NULL`,
// and sets `i` to the matching short option character if possible.
// Otherwise `i` is set to `-1`.
// Note that this procedure uses global state and
// long options only work if `getopt_long` is available.
__attribute__ ((__nonnull__ (2, 3)))
int opt_parse(int, char* const*, char const*, char const* const*);

// The call `opt_parse_str(ptr, str, n)` extracts a string
// from the previously parsed command-line argument,
// writes its index in the array `str` of length `n` into `ptr` and
// returns `true` if possible.
// Otherwise `false` is returned and no copying is done.
__attribute__ ((__nonnull__))
bool opt_parse_str(size_t*, char const* const*, size_t);

// The call `opt_parse_size(ptr, min, max)` extracts a size number
// from the previously parsed command-line argument,
// writes it into `ptr` and returns `true`
// if it is in the closed interval from `min` to `max`.
// Otherwise `false` is returned and no copying is done.
__attribute__ ((__nonnull__))
bool opt_parse_size(size_t*, size_t, size_t);

// The call `opt_parse_fp(ptr, min, max)` extracts a floating-point number
// from the previously parsed command-line argument,
// writes it into `ptr` and returns `true`
// if it is in the closed interval from `min` to `max`.
// Otherwise `false` is returned and no copying is done.
__attribute__ ((__nonnull__))
bool opt_parse_fp(double*, double, double);

#endif
