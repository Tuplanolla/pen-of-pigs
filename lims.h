#ifndef LIMS_H
#define LIMS_H

#include <stddef.h>

// These preprocessor directives define
// the amount of memory reserved for the simulation.
#define DIM_MAX ((size_t) 1 << 2) // 4
#define POLY_MAX ((size_t) 1 << 5) // 32
#define BEAD_MAX ((size_t) 1 << 8) // 256
#define SUBDIV_MAX ((size_t) 1 << 11) // 2.048e+3
#define DIV_MAX ((size_t) 1 << 14) // 16.384e+3
#define STEP_MAX ((size_t) 1 << 32) // 4.294967296e+9

#endif
