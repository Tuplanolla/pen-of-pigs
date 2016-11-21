#include "secs.h"
#include <math.h>
#include <sys/time.h>
#include <time.h>

void secs_to_timeval(struct timeval* const tp, double const s) {
  long int sec = (long int) s;

  tp->tv_sec = sec;
  tp->tv_usec = (long int) ((s - (double) sec) * 1e+6);
}

double secs_from_timeval(struct timeval const* const tp) {
  return (double) tp->tv_sec + (double) tp->tv_usec / 1e+6;
}

void secs_to_timespec(struct timespec* const tp, double const s) {
  long int sec = (long int) s;

  tp->tv_sec = (time_t) sec;
  tp->tv_nsec = (long int) ((s - (double) sec) * 1e+9);
}

double secs_from_timespec(struct timespec const* const tp) {
  return (double) tp->tv_sec + (double) tp->tv_nsec / 1e+9;
}

double secs_now(void) {
  struct timespec tp;
  return clock_gettime(CLOCK_MONOTONIC, &tp) == -1 ?
    (double) NAN : secs_from_timespec(&tp);
}
