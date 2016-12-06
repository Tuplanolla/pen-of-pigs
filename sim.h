#ifndef SIM_H
#define SIM_H

#include "exts.h"
#include <stddef.h>

struct bead;
struct ensem;

// Minimal norm squared from origin.
__attribute__ ((__nonnull__, __pure__))
double bead_norm2(struct ensem const*, struct bead const*);

// Minimal norm from origin.
__attribute__ ((__nonnull__, __pure__))
double bead_norm(struct ensem const*, struct bead const*);

// Minimal distance squared between two points.
__attribute__ ((__nonnull__, __pure__))
double bead_dist2(struct ensem const*, struct bead const*, struct bead const*);

// Minimal distance between two points.
__attribute__ ((__nonnull__, __pure__))
double bead_dist(struct ensem const*, struct bead const*, struct bead const*);

// Magic.
__attribute__ ((__nonnull__))
void sim_run(size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t,
    bool, double, double,
    double (*)(struct ensem const*, struct bead const*, struct bead const*),
    double (*)(struct ensem const*, struct bead const*, struct bead const*),
    double (*)(struct ensem const*, struct bead const*));

#endif
