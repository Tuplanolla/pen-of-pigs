#ifndef SIM_H
#define SIM_H

#include "exts.h"
#include <stddef.h>
#include <stdio.h>

struct bead;
struct ensem;
struct napkin;

// The call `pot_zero(e, r0, r1)` always returns zero.
__attribute__ ((__const__, __nonnull__, __pure__))
double pot_zero(struct ensem const*,
    struct bead const*, struct bead const*);

// The call `potext_zero(e, r)` always returns zero.
__attribute__ ((__const__, __nonnull__, __pure__))
double potext_zero(struct ensem const*, struct bead const*);

// The call `bead_norm2(e, r)` returns the norm squared
// of `r` according to the minimum image convention.
// This is equivalent to `bead_dist2(e, z, r)`, where `z` is the origin.
__attribute__ ((__nonnull__, __pure__))
double bead_norm2(struct ensem const*, struct bead const*);

// The call `bead_norm(e, r)` returns the norm
// of `r` according to the minimum image convention.
// This is equivalent to `sqrt(bead_norm2(e, r))`.
__attribute__ ((__nonnull__, __pure__))
double bead_norm(struct ensem const*, struct bead const*);

// The call `bead_dist2(e, r0, r1)` returns the distance squared
// between `r0` and `r1` according to the minimum image convention.
__attribute__ ((__nonnull__, __pure__))
double bead_dist2(struct ensem const*, struct bead const*, struct bead const*);

// The call `bead_dist(e, r0, r1)` returns the distance
// between `r0` and `r1` according to the minimum image convention.
// This is equivalent to `sqrt(bead_dist2(e, r0, r1))`.
__attribute__ ((__nonnull__, __pure__))
double bead_dist(struct ensem const*, struct bead const*, struct bead const*);

// The call `sim_run(ndim, npoly, nbead, nsubdiv,
// nthrm, nprod, nthrmrec, nprodrec,
// periodic, L, beta, Vint, Vend, Vext)` is magic.
__attribute__ ((__nonnull__))
bool sim_run(size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t,
    bool, double, double,
    double (*)(struct ensem const*, struct bead const*, struct bead const*),
    double (*)(struct ensem const*, struct bead const*, struct bead const*),
    double (*)(struct ensem const*, struct bead const*));

// TODO Organize this mess.

__attribute__ ((__nonnull__, __pure__))
static double K_polybead_bw(struct ensem const*, size_t, size_t);
__attribute__ ((__nonnull__, __pure__))
static double K_polybead_fw(struct ensem const*, size_t, size_t);
__attribute__ ((__nonnull__, __pure__))
static double K_polybead(struct ensem const*, size_t, size_t);
__attribute__ ((__nonnull__, __pure__))
static double K_poly(struct ensem const*, size_t);
__attribute__ ((__nonnull__, __pure__))
static double K_total(struct ensem const*);
__attribute__ ((__nonnull__, __pure__))
static double Vint_bead(struct ensem const*, size_t);
__attribute__ ((__nonnull__, __pure__))
static double Vend_bead(struct ensem const*, size_t);
__attribute__ ((__nonnull__, __pure__))
static double Vext_polybead(struct ensem const*, size_t, size_t);
__attribute__ ((__nonnull__, __pure__))
static double Vext_bead(struct ensem const*, size_t);
__attribute__ ((__nonnull__, __pure__))
static double V_bead(struct ensem const*, size_t);
__attribute__ ((__nonnull__, __pure__))
static double V_total(struct ensem const*);

__attribute__ ((__nonnull__, __pure__))
static double est_pimc_tde(struct ensem const*);
__attribute__ ((__nonnull__, __pure__))
static double est_pigs_crap(struct ensem const*);
__attribute__ ((__nonnull__, __pure__))
static double est_pigs_tde(struct ensem const*);

__attribute__ ((__nonnull__))
static void Rm_const(struct napkin*, double);
__attribute__ ((__nonnull__))
static void Rend_close(struct napkin*);
__attribute__ ((__nonnull__))
static void Rend_open(struct napkin*);
__attribute__ ((__nonnull__))
static void Rend_cycle(struct napkin*);
__attribute__ ((__nonnull__))
static void Rr_rand(struct napkin*);
__attribute__ ((__nonnull__))
static void Rr_randpt(struct napkin*);
__attribute__ ((__nonnull__))
static void Rr_randlatt(struct napkin*);
__attribute__ ((__nonnull__))
static void Rr_ptlatt(struct napkin*);
__attribute__ ((__nonnull__))
static void Rr_circlatt(struct napkin*);

__attribute__ ((__nonnull__))
static void move_accept_ss(struct napkin*);
__attribute__ ((__nonnull__))
static void move_reject_ss(struct napkin*);
__attribute__ ((__nonnull__))
static void move_adjust_ss(struct napkin*);
__attribute__ ((__nonnull__))
static void move_ss(struct napkin*, size_t, size_t);
__attribute__ ((__nonnull__))
static void move_accept_cmd(struct napkin*);
__attribute__ ((__nonnull__))
static void move_reject_cmd(struct napkin*);
__attribute__ ((__nonnull__))
static void move_adjust_cmd(struct napkin*);
__attribute__ ((__nonnull__))
static void move_cmd(struct napkin*, size_t);

__attribute__ ((__noreturn__))
static void move_accept_bisect(struct napkin*);
__attribute__ ((__noreturn__))
static void move_reject_bisect(struct napkin*);
__attribute__ ((__noreturn__))
static void move_bisect(struct napkin*, size_t, size_t);
__attribute__ ((__noreturn__))
static void move_accept_swap(struct napkin*);
__attribute__ ((__noreturn__))
static void move_reject_swap(struct napkin*);
__attribute__ ((__noreturn__))
static void move_swap(struct napkin*, size_t, size_t);

__attribute__ ((__nonnull__))
static double work_ss(struct napkin*);
__attribute__ ((__nonnull__))
static double work_cmd(struct napkin*);

__attribute__ ((__nonnull__))
static void choose(struct napkin*, double);

__attribute__ ((__nonnull__))
static void posdist_accum(struct napkin const*);
__attribute__ ((__nonnull__))
static void paircorr_accum(struct napkin const*);

__attribute__ ((__nonnull__))
static bool res_close(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static FILE* res_open(struct napkin const*, char const*);
__attribute__ ((__nonnull__))
static bool res_use(struct napkin const*, char const*,
    bool (*)(struct napkin const*, FILE*));

__attribute__ ((__nonnull__))
static bool disp_ndim(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_npoly(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_nbead(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_nsubdiv(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_length(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_R_r(struct napkin const*, FILE*, size_t, size_t, size_t);
__attribute__ ((__nonnull__))
static bool disp_pots(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_energy(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_params(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_polys(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_posdist(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_paircorr(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_progress(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_results(struct napkin const*, FILE*);
__attribute__ ((__nonnull__))
static bool disp_wrong_results_fast(struct napkin const*, FILE*);

__attribute__ ((__nonnull__))
static bool prepare_save(struct napkin const*);

__attribute__ ((__nonnull__))
static bool save_const(struct napkin const*);
__attribute__ ((__nonnull__))
static bool save_mut(struct napkin const*);

static void napkin_free(struct napkin*);

__attribute__ ((__malloc__))
static struct napkin* napkin_alloc(size_t, size_t, size_t, size_t,
    size_t, size_t, size_t, size_t);

#endif
