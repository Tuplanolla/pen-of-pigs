#include "apxtime.h"
#include "err.h"
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

static double t0 = 0;

void reset(void) {
  t0 = now();
}

/*
The handling of `errno` respects the following statement from N1570.

> The value of errno may be set to nonzero
> by a library function call whether or not there is an error,
> provided the use of errno is not documented
> in the description of the function in this International Standard.
*/
void warn_with(char const* const func, char const* const file,
    size_t const line, char const* const s) {
  int const n = errno;

  double const t1 = now();

  if (fprintf(stderr, "[%f] (%zu) <%s> %s:%zu: ",
      t1 - t0, (size_t) getpid(), func, file, line) != EOF) {
    errno = n;

    perror(s);
  }
}
