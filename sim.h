#ifndef SIM_H
#define SIM_H

// This file reserves the `bead`, `poly` and `ens` namespaces
// in addition to the `sim` namespace.

#include "exts.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

// This structure represents a bead and is a part of `struct poly`.
struct bead;

// This structure represents a polymer and is a part of `struct ens`.
struct poly;

// This structure represents an ensemble and is a part of `struct sim`.
struct ens;

// This structure contains all the necessary bookkeeping beside `struct ens`.
struct sim;

typedef double (* ens_pot)(struct ens const*,
    struct bead const*, struct bead const*);

// The call `ens_pot_zero(ens, r0, r1)` always returns zero.
__attribute__ ((__const__, __nonnull__, __pure__))
double ens_pot_zero(struct ens const*,
    struct bead const*, struct bead const*);

typedef double (* ens_potext)(struct ens const*, struct bead const*);

// The call `ens_potext_zero(ens, r)` always returns zero.
__attribute__ ((__const__, __nonnull__, __pure__))
double ens_potext_zero(struct ens const*, struct bead const*);

// The call `ens_norm2(ens, r)` returns the norm squared
// of `r` according to the minimum image convention.
// This is equivalent to `ens_dist2(ens, z, r)`, where `z` is the origin.
__attribute__ ((__nonnull__, __pure__))
double ens_norm2(struct ens const*, struct bead const*);

// The call `ens_norm(ens, r)` returns the norm
// of `r` according to the minimum image convention.
// This is equivalent to `sqrt(ens_norm2(ens, r))`.
__attribute__ ((__nonnull__, __pure__))
double ens_norm(struct ens const*, struct bead const*);

// The call `ens_dist2(ens, r0, r1)` returns the distance squared
// between `r0` and `r1` according to the minimum image convention.
__attribute__ ((__nonnull__, __pure__))
double ens_dist2(struct ens const*, struct bead const*, struct bead const*);

// The call `ens_dist(ens, r0, r1)` returns the distance
// between `r0` and `r1` according to the minimum image convention.
// This is equivalent to `sqrt(ens_dist2(ens, r0, r1))`.
__attribute__ ((__nonnull__, __pure__))
double ens_dist(struct ens const*, struct bead const*, struct bead const*);

// TODO Document these kinetic energy calculators.

// Kinetic energy for the degrees of freedom by the equipartition theorem.
__attribute__ ((__nonnull__, __pure__))
double ens_kindf(struct ens const*);

// Kinetic energy within a polymer from one bead to the previous one.
__attribute__ ((__nonnull__, __pure__))
double ens_kin_polybead_bw(struct ens const*, size_t, size_t);

// Kinetic energy within a polymer from one bead to the next one.
// $\lambda = \hbar / 2 m$
// $K = x^2 / 4 \lambda \tau^2 = M x^2$
// $M = m / 2 \hbar \tau^2$
__attribute__ ((__nonnull__, __pure__))
double ens_kin_polybead_fw(struct ens const*, size_t, size_t);

// Kinetic energy within a polymer from each bead to the previous one.
__attribute__ ((__nonnull__, __pure__))
double ens_kin_bead_bw(struct ens const*, size_t);

// Kinetic energy within a polymer from each bead to the next one.
__attribute__ ((__nonnull__, __pure__))
double ens_kin_bead_fw(struct ens const*, size_t);

// Kinetic energy within a polymer from one bead to its two neighbors.
__attribute__ ((__nonnull__, __pure__))
double ens_kin_polybead(struct ens const*, size_t, size_t);

// Kinetic energy within a polymer.
__attribute__ ((__nonnull__, __pure__))
double ens_kin_poly(struct ens const*, size_t);

// Kinetic energy.
__attribute__ ((__nonnull__, __pure__))
double ens_kin_total(struct ens const*);

// TODO Document these kinetic energy calculators.

// Potential energy between certain beads of all polymers.
__attribute__ ((__nonnull__, __pure__))
double ens_potint_bead(struct ens const*, size_t);

// External potential energy for one bead in a polymer.
__attribute__ ((__nonnull__, __pure__))
double ens_potext_polybead(struct ens const*, size_t, size_t);

// External potential energy for every bead in a polymer.
__attribute__ ((__nonnull__, __pure__))
double ens_potext_bead(struct ens const*, size_t);

// Potential energy for one bead over all polymers.
__attribute__ ((__nonnull__, __pure__))
double ens_pot_bead(struct ens const*, size_t);

// Potential energy.
__attribute__ ((__nonnull__, __pure__))
double ens_pot_total(struct ens const*);

// TODO Document these estimators.

typedef double (* sim_est)(struct sim const*, void const*);

// Thermodynamic energy estimator for closed polymers.
// $\langle E_T\rangle = \frac 1{M \tau} \sum_{k = 1}^M
// \Bigl\langle\frac{d N} 2
// - \frac{|R_{k + 1 \bmod M} - R_k|^2}{4 \lambda \tau}
// + \tau V(R_k)\Bigr\rangle$
__attribute__ ((__nonnull__ (1)))
double sim_est_pimc_thermal(struct sim const*, void const*);

// TODO Thermodynamic energy estimator for open polymers.
__attribute__ ((__nonnull__ (1)))
double sim_est_pigs_virial(struct sim const*, void const*);

// TODO Crap energy estimator for open polymers.
__attribute__ ((__nonnull__ (1)))
double sim_est_pigs_mixed(struct sim const*, void const*);

// The call `sim_set_potint(sim, f)` sets `f`
// as the internal (polymer-to-polymer) potential of the simulation `sim`.
__attribute__ ((__nonnull__))
void sim_set_potint(struct sim*, ens_pot);

// The call `sim_set_potext(sim, f)` sets `f`
// as the external (bead-to-origin) potential of the simulation `sim`.
__attribute__ ((__nonnull__))
void sim_set_potext(struct sim*, ens_potext);

// TODO Document these weight configuration setters.

typedef void (* sim_weighter)(struct sim*, void const*);

__attribute__ ((__nonnull__ (1)))
void sim_weight_const(struct sim*, void const*);

// TODO Document these permutation configuration setters.

typedef void (* sim_permer)(struct sim*, void const*);

__attribute__ ((__nonnull__ (1)))
void sim_perm_close(struct sim*, void const*);

__attribute__ ((__nonnull__ (1)))
void sim_perm_open(struct sim*, void const*);

__attribute__ ((__nonnull__ (1)))
void sim_perm_chain(struct sim*, void const*);

__attribute__ ((__nonnull__ (1)))
void sim_perm_random(struct sim*, void const*);

// TODO Document these position configuration setters.

typedef void (* sim_placer)(struct sim*, size_t, struct bead const*, double,
    void const*);

__attribute__ ((__nonnull__ (1, 3)))
void sim_placer_point(struct sim*, size_t, struct bead const*, double,
    void const*);

__attribute__ ((__nonnull__ (1, 3)))
void sim_placer_random(struct sim*, size_t, struct bead const*, double,
    void const*);

__attribute__ ((__nonnull__ (1, 3)))
void sim_placer_knot(struct sim*, size_t, struct bead const*, double,
    void const*);

__attribute__ ((__nonnull__ (1, 2)))
void sim_place_point(struct sim*, sim_placer, void const*);

__attribute__ ((__nonnull__ (1, 2)))
void sim_place_random(struct sim*, sim_placer, void const*);

__attribute__ ((__nonnull__ (1, 2)))
void sim_place_lattice(struct sim*, sim_placer, void const*);

// TODO Document these movers.

typedef void (* sim_mover)(struct sim*);

__attribute__ ((__nonnull__))
void sim_move_null(struct sim*);

__attribute__ ((__nonnull__))
void sim_move_accept_ss(struct sim*);

__attribute__ ((__nonnull__))
void sim_move_reject_ss(struct sim*);

__attribute__ ((__nonnull__))
void sim_move_adjust_ss(struct sim*);

__attribute__ ((__nonnull__))
void sim_move_ss(struct sim*, size_t, size_t);

__attribute__ ((__nonnull__))
void sim_move_accept_cmd(struct sim*);

__attribute__ ((__nonnull__))
void sim_move_reject_cmd(struct sim*);

__attribute__ ((__nonnull__))
void sim_move_adjust_cmd(struct sim*);

__attribute__ ((__nonnull__))
void sim_move_cmd(struct sim*, size_t);

__attribute__ ((__noreturn__))
void sim_move_accept_bisect(struct sim*);

__attribute__ ((__noreturn__))
void sim_move_reject_bisect(struct sim*);

__attribute__ ((__noreturn__))
void sim_move_bisect(struct sim*, size_t, size_t);

__attribute__ ((__noreturn__))
void sim_move_accept_swap(struct sim*);

__attribute__ ((__noreturn__))
void sim_move_reject_swap(struct sim*);

__attribute__ ((__noreturn__))
void sim_move_swap(struct sim*, size_t, size_t);

// TODO Document these move proposers.

typedef double (* sim_proposer)(struct sim*);

__attribute__ ((__nonnull__))
double sim_propose_ss(struct sim*);

__attribute__ ((__nonnull__))
double sim_propose_cmd(struct sim*);

// TODO Document these move deciders.

typedef void (* sim_decider)(struct sim*);

__attribute__ ((__nonnull__))
void sim_decide_mq(struct sim*, double);

// TODO Document these resource managers.

typedef bool (* sim_printer)(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__))
bool sim_res_close(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
FILE* sim_res_open(struct sim const*, char const*);

__attribute__ ((__nonnull__ (1, 2, 3)))
bool sim_res_print(struct sim const*, char const*, sim_printer, void const*);

// TODO Document these resource printers.

__attribute__ ((__nonnull__ (1, 2)))
bool print_length(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_periodic(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_pots(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_ndim(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_npoly(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_nbead(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_nsubdiv(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_ndiv(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_energy_bead(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_energy(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_energy_corrtime(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_params(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_posdist(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_raddist(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_polys(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_progress(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_results(struct sim const*, FILE*, void const*);

__attribute__ ((__nonnull__ (1, 2)))
bool print_wrong_results_fast(struct sim const*, FILE*, void const*);

// TODO Document these environment managers.

__attribute__ ((__nonnull__))
bool sim_init_fs(struct sim const*);

__attribute__ ((__nonnull__))
bool sim_fini_fs(struct sim const*);

__attribute__ ((__nonnull__))
bool sim_save_const(struct sim const*);

__attribute__ ((__nonnull__))
bool sim_save_mut(struct sim const*);

// TODO Document these simulation managers.

void sim_free(struct sim*);

__attribute__ ((__malloc__))
struct sim* sim_alloc(size_t, size_t, size_t, size_t,
    size_t, size_t, size_t, size_t, size_t,
    bool, bool, double, double, double);

// The call `sim_run(sim)` is magic.
__attribute__ ((__nonnull__))
bool sim_run(struct sim*);

#endif
