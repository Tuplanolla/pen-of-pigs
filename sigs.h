#ifndef SIGS_H
#define SIGS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define SIGS_MIN (SIG_ATOMIC_MIN + 2)
#define SIGS_MAX (SIG_ATOMIC_MAX - 1)

// The call `sigs_handler(signum)` handles the signal with the number `signum`.
// This can and should be registered as a signal handler with `sigs_register`.
void sigs_handler(int);

// The statement `i = sigs_register(sigs, n)` tries to register `sigs_handler`
// as the signal handler of each signal in the array `sigs` of length `n`.
// If the operation is successful, `i` is set to `SIZE_MAX`.
// Otherwise the registering is stopped at the first failure and
// `i` is set the index of the first signal that failed to be registered.
size_t sigs_register(int const*, size_t);

// This procedure returns `true` if a signal has not been caught.
// Otherwise `false` is returned.
bool sigs_unset(void);

// This procedure returns `true` if a signal has been caught.
// Otherwise `false` is returned.
bool sigs_set(void);

// This procedure returns `true`
// if a signal below the representable minimum (negative) has been caught.
// Otherwise `false` is returned.
bool sigs_under(void);

// The call `sigs_normal(ptr)` returns `true`
// if a normal signal (with a small enough magnitude) has been caught and
// copies its value into `ptr` if `ptr` is not `NULL`.
// Otherwise `false` is returned and no copying is done.
bool sigs_normal(int*);

// This procedure returns `true` if a null signal (zero) has been caught.
// Otherwise `false` is returned.
bool sigs_null(void);

// This procedure returns `true`
// if a signal above the representable maximum (positive) has been caught.
// Otherwise `false` is returned.
bool sigs_over(void);

// The call `sigs_forget()` forgets about any previously caught signals.
// Calling `sigs_forget` before using any other signal handling procedures
// from this compilation unit is not necessary.
void sigs_forget(void);

// The statement `p = sigs_use(ptr)` is equivalent
// to `p = sigs_normal(ptr)` followed by `sigs_forget()` atomically.
bool sigs_use(int*);

#endif
