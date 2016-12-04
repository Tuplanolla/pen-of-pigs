#ifndef EXTS_H
#define EXTS_H

#include <stdbool.h>

// This disables GNU extensions if they are unsupported.
#if !defined __GNUC__ || __GNUC__ < 4
#ifndef __attribute__
#define __attribute__(_)
#endif
#endif

// The preprocessor directives `BEGIN` and `END` allow converting
// multiple statements into a single statement
// with local scope and no redundant semicolons.
// They must always appear together.
#define BEGIN do {
#define END } while (false)

// This preprocessor directive imitates `static_assert` if it is not available.
// Due to technical limitations,
// each `static_assert` must be on its own line to avoid naming conflicts.
#ifndef static_assert
#define _static_assert_line(p, n) \
  __attribute__ ((__unused__)) \
  static int const _static_assert_##n[(p) ? 1 : -1]
#define _static_assert(p, n) _static_assert_line((p), n)
#define static_assert(p, _) _static_assert((p), __LINE__)
#endif

// The preprocessor directive `dynamic_assert(p, s)`
// is equivalent to `assert(p)`,
// which makes it consistent with `static_assert`.
#ifndef dynamic_assert
#define dynamic_assert(p, s) assert(p)
#endif

// This ensures that exactly one of `NDEBUG` or `DEBUG` is defined.
#ifdef DEBUG
#ifdef NDEBUG
static_assert(false, "contradictory debug directives");
#else
#endif
#else
#ifndef NDEBUG
#define NDEBUG
#endif
#endif

// This preprocessor directive makes it possible to write `for ever`.
// This is very important.
#ifndef ever
#define ever (;;)
#endif

#endif
