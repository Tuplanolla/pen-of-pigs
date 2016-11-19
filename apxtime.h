#ifndef APXTIME_H
#define APXTIME_H

#include "exts.h"
#include <sys/time.h>
#include <time.h>

/*
The call `pack_timeval(tp, s)` sets
the time structure `tp` to approximately `s` seconds.
*/
__attribute__ ((__nonnull__))
void pack_timeval(struct timeval*, double);

/*
The call `unpack_timeval(tp)` returns
the approximate number of seconds in the time structure `tp`.
*/
__attribute__ ((__nonnull__))
double unpack_timeval(struct timeval const* const);

/*
The call `pack_timespec(tp, s)` sets
the time structure `tp` to approximately `s` seconds.
*/
__attribute__ ((__nonnull__))
void pack_timespec(struct timespec*, double);

/*
The call `unpack_timespec(tp)` returns
the approximate number of seconds in the time structure `tp`.
*/
__attribute__ ((__nonnull__))
double unpack_timespec(struct timespec const* const);

/*
The call `now()` returns
the approximate time right now or `NAN` on error.
*/
double now(void);

#endif
