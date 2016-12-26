#include "sigs.h"
#include <limits.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>

#define SIGS_UNSET (SIGS_MIN - 2)
#define SIGS_UNDER (SIGS_MIN - 1)
#define SIGS_NULL 0
#define SIGS_OVER (SIGS_MAX + 1)

static volatile sig_atomic_t sig = SIGS_UNSET;

static bool signum_under(int const signum) {
  return signum < SIGS_MIN;
}

static bool signum_normal(int const signum) {
  return signum >= SIGS_MIN && signum <= SIGS_MAX;
}

static bool signum_null(int const signum) {
  return signum == 0;
}

static bool signum_over(int const signum) {
  return signum > SIGS_MAX;
}

void sigs_handler(int const signum) {
  sig = signum_normal(signum) ? (sig_atomic_t) signum :
    signum_null(signum) ? SIGS_NULL :
    signum_over(signum) ? SIGS_OVER :
    signum_under(signum) ? SIGS_UNDER :
    SIGS_UNSET;
}

// This is possible with `signal`
// if it provides BSD semantics instead of System V semantics,
// but why bother with `signal` when `sigaction` exists?
size_t sigs_register(int const* const sigs, size_t const n) {
  struct sigaction act;
  act.sa_handler = sigs_handler;
  sigfillset(&act.sa_mask);
  act.sa_flags = SA_RESTART;

  for (size_t i = 0; i < n; ++i)
    if (sigaction(sigs[i], &act, NULL) == -1)
      return i;

  return SIZE_MAX;
}

size_t sigs_unregister(int const* const sigs, size_t const n) {
  struct sigaction act;
  act.sa_handler = SIG_DFL;
  sigfillset(&act.sa_mask);
  act.sa_flags = 0;

  for (size_t i = 0; i < n; ++i)
    if (sigaction(sigs[i], &act, NULL) == -1)
      return i;

  return SIZE_MAX;
}

bool sigs_unset(void) {
  return sig == SIGS_UNSET;
}

bool sigs_set(void) {
  return sig != SIGS_UNSET;
}

bool sigs_under(void) {
  return sig == SIGS_MIN;
}

bool sigs_normal(int* const ptr) {
  sig_atomic_t const n = sig;

  bool const p = n >= SIGS_MIN && n <= SIGS_MAX;

  if (p && ptr != NULL)
    *ptr = (int) n;

  return p;
}

bool sigs_null(void) {
  return sig == SIGS_NULL;
}

bool sigs_over(void) {
  return sig == SIGS_MAX;
}

void sigs_forget(void) {
  sig = SIGS_UNSET;
}

bool sigs_use(int* const ptr) {
  sigset_t set;
  sigfillset(&set);

  sigset_t oldset;
  bool const p = sigprocmask(SIG_SETMASK, &set, &oldset) == -1;

  bool const q = sigs_normal(ptr);
  if (q)
    sigs_forget();

  if (!p)
    (void) sigprocmask(SIG_SETMASK, &oldset, NULL);

  return q;
}
