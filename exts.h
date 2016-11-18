#ifndef EXTS_H
#define EXTS_H

/*
This disables GNU extensions if they are unsupported.
*/
#if !defined __GNUC__ || __GNUC__ < 4
#ifndef __attribute__
#define __attribute__(_)
#endif
#endif

/*
This preprocessor directive imitates `static_assert` if it is not available.
Note that each `static_assert` must be on its own line to avoid conflicts.
*/
#ifndef static_assert
#define _static_assert_line(p, n) \
  __attribute__ ((__unused__)) \
  static int const _static_assert_##n[(p) ? 1 : -1]
#define _static_assert(p, n) _static_assert_line((p), n)
#define static_assert(p, _) _static_assert((p), __LINE__)
#endif

/*
This preprocessor directive allows one to write `for ever`.
*/
#ifndef ever
#define ever (;;)
#endif

#endif
