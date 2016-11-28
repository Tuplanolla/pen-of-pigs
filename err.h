#ifndef ERR_H
#define ERR_H

#include "exts.h"
#include <stdbool.h>
#include <stdlib.h>

// The call `err_reset()` initializes the timing mechanism
// `err_msg` and `err_abort` are based on.
// Calling `err_reset` before using any other error handling procedures
// from this compilation unit is not necessary.
// If `err_reset` is never called,
// the reference point for timestamps in error messages is chosen arbitrarily.
void err_reset(void);

// The call `err_msg_with(func, file, line, str)` prints a detailed error message
// from the procedure `func` on line `line` of file `file`
// to the standard error output stream,
// where the pointer `str` may be `NULL` or point to a string
// that contains the name of the procedure or variable that caused the error.
// The error number is taken from the external variable `errno`,
// which may be reset in the process.
// The format of the error message consists of
//
// * `"[time]"` and `"(pid)"`, mimicking `dmesg`;
// * `"<procedure>"`, mimicking `objdump`;
// * `"file:line: message"`, mimicking `cc` and
// * `"procedure: message"`, mimicking `perror`.
//
// It is better to use `err_msg` instead of calling `err_msg_with` directly.
__attribute__ ((__nonnull__ (1, 2)))
void err_msg_with(char const*, char const*, size_t, char const*);

// The preprocessor directive `err_msg(p)` prints a detailed error message
// to the standard error output stream,
// where the pointer `p` may be `NULL` or
// point to the procedure or variable that caused the error.
// See `err_msg_with(p)` for details.
#define err_msg(p) BEGIN \
  err_msg_with(__func__, __FILE__, (size_t) __LINE__, (p) == NULL ? NULL : #p); \
END

// The call `err_abort(p)` is equivalent to `err_msg(p)` followed by `abort()`.
#define err_abort(p) BEGIN \
  err_msg(p); \
  abort(); \
END

#endif
