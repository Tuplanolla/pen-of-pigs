#include "sig.h"
#include <limits.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>

#define LOCK_OFF 0
#define LOCK_ON 1

#define SIG_UNSET (SIG_MIN - 2)
#define SIG_UNDER (SIG_MIN - 1)
#define SIG_NULL 0
#define SIG_OVER (SIG_MAX + 1)

static volatile sig_atomic_t lock = LOCK_OFF;
static volatile sig_atomic_t sig = SIG_UNSET;

static bool signum_under(int const signum) {
  return signum < SIG_MIN;
}

static bool signum_normal(int const signum) {
  return signum >= SIG_MIN && signum <= SIG_MAX;
}

static bool signum_null(int const signum) {
  return signum == 0;
}

static bool signum_over(int const signum) {
  return signum > SIG_MAX;
}

void sig_handler(int const signum) {
  while (lock == LOCK_ON);

  sig = signum_normal(signum) ? (sig_atomic_t) signum :
    signum_null(signum) ? SIG_NULL :
    signum_over(signum) ? SIG_OVER :
    signum_under(signum) ? SIG_UNDER :
    SIG_UNSET;
}

size_t sig_register(int const* const xs, size_t const n) {
  struct sigaction sa;
  sa.sa_handler = sig_handler;
  sigfillset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART;

  for (size_t i = 0;
      i < n;
      ++i)
    if (sigaction(xs[i], &sa, NULL) == -1)
      return i;

  return SIZE_MAX;
}

bool sig_unset(void) {
  return sig == SIG_UNSET;
}

bool sig_set(void) {
  return sig != SIG_UNSET;
}

bool sig_under(void) {
  return sig == SIG_MIN;
}

bool sig_normal(int* const signum) {
  sig_atomic_t const n = sig;

  bool const p = n >= SIG_MIN && n <= SIG_MAX;

  if (p && signum != NULL)
    *signum = (int) n;

  return p;
}

bool sig_null(void) {
  return sig == SIG_NULL;
}

bool sig_over(void) {
  return sig == SIG_MAX;
}

void sig_forget(void) {
  sig = SIG_UNSET;
}

bool sig_use(int* const signum) {
  lock = LOCK_ON;

  bool const p = sig_normal(signum);
  if (p)
    sig_forget();

  lock = LOCK_OFF;

  return p;
}
