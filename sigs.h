#ifndef SIGS_H
#define SIGS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define SIGS_MIN (SIG_ATOMIC_MIN + 2)
#define SIGS_MAX (SIG_ATOMIC_MAX - 1)

/*
The call `sigs_handler(signum)` handles the signal with the number `signum`.
This can and should be registered as a signal handler with `sigs_register`.
*/
void sigs_handler(int);

/*
The call `i = sigs_register(xs, n)` registers `sigs_handler`
for each signal in the array `xs` of length `n`.
If the operation is successful, `i` becomes `SIZE_MAX`.
Otherwise `i` becomes the index of the first signal that failed to register and
registering the rest of the signals is not attempted.
*/
size_t sigs_register(int const*, size_t);

/*
This procedure returns `true` if a signal has not been caught.
*/
bool sigs_unset(void);

/*
This procedure returns `true` if a signal has been caught.
*/
bool sigs_set(void);

/*
This procedure returns `true` if
a signal below the representable minimum (negative) has been caught.
*/
bool sigs_under(void);

/*
The call `sigs_normal(signum)` returns `true` if
a normal signal (with a small magnitude) has been caught and
copies its value into `signum` if it is not `NULL`.
Otherwise `false` is returned.
*/
bool sigs_normal(int*);

/*
This procedure returns `true` if a null signal (zero) has been caught.
*/
bool sigs_null(void);

/*
This procedure returns `true` if
a signal above the representable maximum (positive) has been caught.
*/
bool sigs_over(void);

/*
This procedure forgets any previously caught signals.
Calling `sigs_forget` before the first use is not necessary.
*/
void sigs_forget(void);

/*
This procedure calls `sigs_normal` and then `sigs_forget`.
*/
bool sigs_use(int*);

#endif
