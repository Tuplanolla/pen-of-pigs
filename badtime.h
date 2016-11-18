#ifndef BADTIME_H
#define BADTIME_H

#include "exts.h"
#include <sys/time.h>
#include <time.h>

/*
The call `pack_timeval(tp, s)` sets
the time structure `tp` to approximately `s` seconds.
*/
__attribute__ ((__nonnull__))
void pack_timeval(struct timeval* tp, double s);

/*
The call `unpack_timeval(tp)` returns
the approximate number of seconds in the time structure `tp`.
*/
__attribute__ ((__nonnull__))
double unpack_timeval(struct timeval const* const tp);

/*
The call `pack_timespec(tp, s)` sets
the time structure `tp` to approximately `s` seconds.
*/
__attribute__ ((__nonnull__))
void pack_timespec(struct timespec* tp, double s);

/*
The call `unpack_timespec(tp)` returns
the approximate number of seconds in the time structure `tp`.
*/
__attribute__ ((__nonnull__))
double unpack_timespec(struct timespec const* const tp);

/*
The call `now()` returns
the approximate time right now or `NAN` on error.
*/
double now(void);

#endif
