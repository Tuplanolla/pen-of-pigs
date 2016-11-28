#ifndef SECS_H
#define SECS_H

#include "exts.h"
#include <sys/time.h>
#include <time.h>

// The call `secs_to_timeval(tp, s)` sets
// the time structure `tp` to approximately `s` seconds.
__attribute__ ((__nonnull__))
void secs_to_timeval(struct timeval*, double);

// The call `secs_from_timeval(tp)` returns
// the approximate number of seconds in the time structure `tp`.
__attribute__ ((__nonnull__))
double secs_from_timeval(struct timeval const*);

// The call `secs_to_timespec(tp, s)` sets
// the time structure `tp` to approximately `s` seconds.
__attribute__ ((__nonnull__))
void secs_to_timespec(struct timespec*, double);

// The call `secs_from_timespec(tp)` returns
// the approximate number of seconds in the time structure `tp`.
__attribute__ ((__nonnull__))
double secs_from_timespec(struct timespec const*);

// The call `secs_now()` returns
// the approximate monotonic time in seconds right now if possible.
// Otherwise `NAN` is returned.
double secs_now(void);

#endif
