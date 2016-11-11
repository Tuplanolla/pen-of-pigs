#ifndef REPORT_H
#define REPORT_H

#include <stdbool.h>
#include <stdlib.h>

/*
The preprocessor directives `BEGIN` and `END` allow converting
multiple statements into a single statement
with local scope and no redundant semicolons.
*/
#define BEGIN do {
#define END } while (false)

/*
The call `reset()` initializes the timing mechanism `warn` is based on.
Calling `reset` is optional; if `reset` is never called,
the reference point for timestamps in error messages is chosen arbitrarily.
*/
void reset(void);

/*
The call `warn_with(func, file, line, s)` prints a detailed error message
from the procedure `func` on line `line` of file `file`, where the pointer `s`
is `NULL` or points to a string with the name of the procedure
that caused the error.
You should use `warn` instead of calling `warn_with` directly.
*/
void warn_with(char const* func, char const* file, size_t line, char const* s);

/*
The call `warn(p)` prints a detailed error message
to the standard error output stream, where the pointer `p`
is `NULL` or points to the procedure that caused the error.
The error number is taken from the external variable `errno`,
which may be reset in the process.

The format of the error message consists of

* `"[time]"` and `"(pid)"`, mimicking `dmesg`;
* `"<procedure>"`, mimicking `objdump`;
* `"file:line: message"`, mimicking `cc` and
* `"procedure: message"`, mimicking `perror`.
*/
#define warn(p) BEGIN \
  warn_with(__func__, __FILE__, (size_t) __LINE__, (p) == NULL ? NULL : #p); \
END

/*
The call `halt(p)` is equivalent to `warn(p)` followed by `abort()`.
*/
#define halt(p) BEGIN \
  warn(p); \
  abort(); \
END

#endif
