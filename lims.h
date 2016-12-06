#ifndef LIMS_H
#define LIMS_H

#include <stddef.h>

// These preprocessor directives define
// the amount of memory reserved for the simulation.
#define DIM_MAX ((size_t) 1 << 2) // 4
#define POLY_MAX ((size_t) 1 << 5) // 32
#define BEAD_MAX ((size_t) 1 << 8) // 256
#define SUBDIV_MAX ((size_t) 1 << 11) // 2.048e+3

#endif
