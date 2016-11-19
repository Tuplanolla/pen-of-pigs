#ifndef SIG_H
#define SIG_H

#include <stdbool.h>
#include <stdint.h>

#define SIG_MIN (SIG_ATOMIC_MIN + 2)
#define SIG_MAX (SIG_ATOMIC_MAX - 1)

/*
The call `sig_handler(signum)` handles the signal with the number `signum`.
This can and should be registered as a signal handler with `sig_register`.
*/
void sig_handler(int);

/*
The call `i = sig_register(xs, n)` registers `sig_handler`
for each signal in the array `xs` of length `n`.
If the operation is successful, `i` is `SIZE_MAX`.
Otherwise `i` becomes the index of the first failure and
the rest are not attempted.
*/
size_t sig_register(int const*, size_t);

/*
This procedure returns `true` if a signal has not been caught.
*/
bool sig_unset(void);

/*
This procedure returns `true` if a signal has been caught.
*/
bool sig_set(void);

/*
This procedure returns `true` if
a signal below the representable minimum (negative) has been caught.
*/
bool sig_under(void);

/*
The call `sig_normal(signum)` returns `true` if
a normal signal (with a small magnitude) has been caught and
copies its value into `signum` if it is not `NULL`.
*/
bool sig_normal(int*);

/*
This procedure returns `true` if a null signal (zero) has been caught.
*/
bool sig_null(void);

/*
This procedure returns `true` if
a signal above the representable maximum (positive) has been caught.
*/
bool sig_over(void);

/*
This procedure forgets any previously caught signals.
Calling `sig_forget` before the first use is not necessary.
*/
void sig_forget(void);

/*
This procedure combines `sig_normal` and `sig_foget`.
*/
bool sig_use(int*);

#endif
