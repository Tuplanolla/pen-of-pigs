#ifndef EXTS_H
#define EXTS_H

/*
This preprocessor directive allows one to write `for ever`.
*/
#define ever (;;)

/*
This disables GNU extensions if they are unsupported.
*/
#if !defined __GNUC__ || __GNUC__ < 4
#define __attribute__(_)
#endif

#endif
