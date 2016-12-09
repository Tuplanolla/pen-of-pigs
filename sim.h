#ifndef SIM_H
#define SIM_H

#include "exts.h"
#include <stddef.h>
#include <stdio.h>

struct bead;
struct ensem;
struct sim;

// The call `sim_pot_zero(ensem, r0, r1)` always returns zero.
__attribute__ ((__const__, __nonnull__, __pure__))
double sim_pot_zero(struct ensem const*,
    struct bead const*, struct bead const*);

// The call `sim_potext_zero(ensem, r)` always returns zero.
__attribute__ ((__const__, __nonnull__, __pure__))
double sim_potext_zero(struct ensem const*, struct bead const*);

// The call `sim_norm2(ensem, r)` returns the norm squared
// of `r` according to the minimum image convention.
// This is equivalent to `sim_dist2(ensem, z, r)`, where `z` is the origin.
__attribute__ ((__nonnull__, __pure__))
double sim_norm2(struct ensem const*, struct bead const*);

// The call `sim_norm(ensem, r)` returns the norm
// of `r` according to the minimum image convention.
// This is equivalent to `sqrt(sim_norm2(ensem, r))`.
__attribute__ ((__nonnull__, __pure__))
double sim_norm(struct ensem const*, struct bead const*);

// The call `sim_dist2(ensem, r0, r1)` returns the distance squared
// between `r0` and `r1` according to the minimum image convention.
__attribute__ ((__nonnull__, __pure__))
double sim_dist2(struct ensem const*, struct bead const*, struct bead const*);

// The call `sim_dist(ensem, r0, r1)` returns the distance
// between `r0` and `r1` according to the minimum image convention.
// This is equivalent to `sqrt(sim_dist2(ensem, r0, r1))`.
__attribute__ ((__nonnull__, __pure__))
double sim_dist(struct ensem const*, struct bead const*, struct bead const*);

// The call `sim_get_ensem(sim)` returns the ensemble of the simulation `sim`.
struct ensem* sim_get_ensem(struct sim*);

// The call `sim_periodic(ens, p)` sets the periodicity
// of the ensemble `ens` to `p`.
__attribute__ ((__nonnull__))
void sim_periodic(struct ensem*, bool);

typedef double (* sim_pot)(struct ensem const*,
    struct bead const*, struct bead const*);

typedef double (* sim_potext)(struct ensem const*, struct bead const*);

// The call `sim_set_potint(ens, f)` sets `f`
// as the internal (polymer-to-polymer) potential of the ensemble `ens`.
__attribute__ ((__nonnull__))
void sim_set_potint(struct ensem*, sim_pot);

// The call `sim_set_potend(ens, f)` sets `f`
// as the additional (end-to-end) potential of the ensemble `ens`.
__attribute__ ((__nonnull__))
void sim_set_potend(struct ensem*, sim_pot);

// The call `sim_set_potext(ens, f)` sets `f`
// as the external (bead-to-origin) potential of the ensemble `ens`.
__attribute__ ((__nonnull__))
void sim_set_potext(struct ensem*, sim_potext);

// TODO Organize this mess.

__attribute__ ((__nonnull__, __pure__))
static double sim_kin_polybead_bw(struct ensem const*, size_t, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_kin_polybead_fw(struct ensem const*, size_t, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_kin_polybead(struct ensem const*, size_t, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_kin_poly(struct ensem const*, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_kin_total(struct ensem const*);

__attribute__ ((__nonnull__, __pure__))
static double sim_potint_bead(struct ensem const*, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_potend_bead(struct ensem const*, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_potext_polybead(struct ensem const*, size_t, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_potext_bead(struct ensem const*, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_pot_bead(struct ensem const*, size_t);

__attribute__ ((__nonnull__, __pure__))
static double sim_pot_total(struct ensem const*);

__attribute__ ((__nonnull__, __pure__))
static double est_pimc_tde(struct ensem const*);

__attribute__ ((__nonnull__, __pure__))
static double est_pigs_crap(struct ensem const*);

__attribute__ ((__nonnull__, __pure__))
static double est_pigs_tde(struct ensem const*);

__attribute__ ((__nonnull__))
void sim_mass_const(struct sim*, double);

__attribute__ ((__nonnull__))
void sim_perm_close(struct sim*);

__attribute__ ((__nonnull__))
void sim_perm_open(struct sim*);

__attribute__ ((__nonnull__))
void sim_perm_chain(struct sim*);

__attribute__ ((__nonnull__))
void sim_perm_random(struct sim*);

typedef void (* sim_placer)(struct sim*,
    size_t, struct bead const*, double);

__attribute__ ((__nonnull__))
void sim_placer_point(struct sim*, size_t, struct bead const*, double);

__attribute__ ((__nonnull__))
void sim_placer_random(struct sim*, size_t, struct bead const*, double);

__attribute__ ((__nonnull__))
void sim_placer_unknot(struct sim*, size_t, struct bead const*, double);

__attribute__ ((__nonnull__))
void sim_placer_tinyuk(struct sim*, size_t, struct bead const*, double);

__attribute__ ((__nonnull__))
void sim_place_point(struct sim*, sim_placer);

__attribute__ ((__nonnull__))
void sim_place_random(struct sim*, sim_placer);

__attribute__ ((__nonnull__))
void sim_place_lattice(struct sim*, sim_placer);

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

typedef double (* sim_proposer)(struct sim*);

__attribute__ ((__nonnull__))
double sim_propose_ss(struct sim*);

__attribute__ ((__nonnull__))
double sim_propose_cmd(struct sim*);

typedef void (* sim_decider)(struct sim*);

__attribute__ ((__nonnull__))
static void decide(struct sim*, double);

__attribute__ ((__nonnull__))
static void posdist_accum(struct sim const*);

__attribute__ ((__nonnull__))
static void paircorr_accum(struct sim const*);

__attribute__ ((__nonnull__))
static void ensem_extents(struct ensem const*, double*, double*);

typedef bool (* sim_printer)(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool res_close(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static FILE* res_open(struct sim const*, char const*);

__attribute__ ((__nonnull__))
static bool res_print(struct sim const*, char const*, sim_printer);

__attribute__ ((__nonnull__))
static bool print_length(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_periodic(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_pots(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_ndim(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_npoly(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_nbead(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_nsubdiv(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_energy(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_energy_corrtime(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_params(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_posdist(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_paircorr(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_polys1(struct sim const*, FILE*, size_t, size_t, size_t);

__attribute__ ((__nonnull__))
static bool print_polys(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_progress(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_results(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool print_wrong_results_fast(struct sim const*, FILE*);

__attribute__ ((__nonnull__))
static bool prepare_save(struct sim const*);

__attribute__ ((__nonnull__))
static bool save_const(struct sim const*);

__attribute__ ((__nonnull__))
static bool save_mut(struct sim const*);

void sim_free(struct sim*);

__attribute__ ((__malloc__))
struct sim* sim_alloc(size_t, size_t, size_t, size_t,
    size_t, size_t, size_t, size_t,
    double, double, double);

// The call `sim_run(sim)` is magic.
__attribute__ ((__nonnull__))
bool sim_run(struct sim*);

#endif
