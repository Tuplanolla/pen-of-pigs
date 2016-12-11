#include "err.h"
#include "exts.h"
#include "fp.h"
#include "lims.h"
#include "ran.h"
#include "sigs.h"
#include "sim.h"
#include "size.h"
#include "stats.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _GNU_SOURCE
#ifdef DEBUG
#include <fenv.h>
#endif
#endif

// This structure contains the numbers of allocated objects by type.
// The numbers never exceed the statically defined limits
// `DIM_MAX`, `POLY_MAX`, `BEAD_MAX` and `SUBDIV_MAX`.
struct memb {
  size_t dim;
  size_t poly;
  size_t bead;
  size_t subdiv;
};

// This structure contains the numbers of total and recorded
// thermalization and production steps.
// The numbers of recording steps never exceed
// that of the corresponding total steps.
struct step {
  size_t thrm;
  size_t prod;
  size_t thrmrec;
  size_t prodrec;
};

// This structure represents a bead and is a part of `struct poly`.
struct bead {
  double d[DIM_MAX];
};

// This structure represents a polymer and is a part of `struct ensem`.
struct poly {
  double m;
  size_t ifrom;
  size_t ito;
  struct bead r[BEAD_MAX];
};

// This structure represents an ensemble and is a part of `struct sim`.
struct ensem {
  bool periodic;
  double L;
  double tau;
  sim_pot Vint;
  sim_pot Vend;
  sim_potext Vext;
  struct memb nmemb;
  struct step nstep;
  struct step istep;
  struct poly R[POLY_MAX];
};

// This structure contains all the necessary bookkeeping.
struct sim {
  gsl_rng* rng;
  struct stats* tde;
  size_t* p;
  size_t* g;
  sim_decider accept;
  sim_decider reject;
  sim_decider adjust;
  struct {
    struct {
      size_t ipoly;
      size_t ibead;
      struct bead r;
    } ssm;
    struct {
      size_t ipoly;
      struct poly R;
    } cmd;
  } hist;
  struct {
    struct {
      size_t acc;
      size_t rej;
      double h;
    } ssm;
    struct {
      size_t acc;
      size_t rej;
      double h;
    } cmd;
  } params;
  struct ensem ens;
  char id[BUFSIZ];
};

double sim_pot_zero(
    __attribute__ ((__unused__)) struct ensem const* const ens,
    __attribute__ ((__unused__)) struct bead const* const r0,
    __attribute__ ((__unused__)) struct bead const* const r1) {
  return 0.0;
}

double sim_potext_zero(
    __attribute__ ((__unused__)) struct ensem const* const ens,
    __attribute__ ((__unused__)) struct bead const* const r) {
  return 0.0;
}

double sim_norm2(struct ensem const* const ens,
    struct bead const* const r) {
  double s = 0.0;

  double (* const f)(double, double) =
    ens->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ens->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r->d[idim], ens->L));

  return s;
}

double sim_norm(struct ensem const* const ens,
    struct bead const* const r) {
  return sqrt(sim_norm2(ens, r));
}

double sim_dist2(struct ensem const* const ens,
    struct bead const* const r0, struct bead const* const r1) {
  double s = 0.0;

  double (* const f)(double, double) =
    ens->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ens->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r1->d[idim] - r0->d[idim], ens->L));

  return s;
}

double sim_dist(struct ensem const* const ens,
    struct bead const* const r0, struct bead const* const r1) {
  return sqrt(sim_dist2(ens, r0, r1));
}

struct ensem* sim_get_ensem(struct sim* const sim) {
  return &sim->ens;
}

void sim_periodic(struct ensem* const ens, bool const p) {
  ens->periodic = p;
}

void sim_set_potint(struct ensem* const ens,
    double (* const V)(struct ensem const*,
      struct bead const*, struct bead const*)) {
  ens->Vint = V;
}

void sim_set_potend(struct ensem* const ens,
    double (* const V)(struct ensem const*,
      struct bead const*, struct bead const*)) {
  ens->Vend = V;
}

void sim_set_potext(struct ensem* const ens,
    double (* const V)(struct ensem const*, struct bead const*)) {
  ens->Vext = V;
}

// Kinetic energy within a polymer from one bead to the previous one.
static double sim_kin_polybead_bw(struct ensem const* const ens,
    size_t const ipoly, size_t const ibead) {
  double const M = ens->R[ipoly].m / (2.0 * gsl_pow_2(ens->tau));

  if (ibead == 0) {
    size_t const jpoly = ens->R[ipoly].ifrom;
    size_t const jbead = ens->nmemb.bead - 1;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return M * sim_dist2(ens,
          &ens->R[ipoly].r[ibead], &ens->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead - 1;

    return M * sim_dist2(ens,
        &ens->R[ipoly].r[ibead], &ens->R[ipoly].r[jbead]);
  }
}

// Kinetic energy within a polymer from one bead to the next one.
// \lambda = \hbar / 2 m, K = x^2 / 4 \lambda \tau^2 = M x^2
// M = m / 2 \hbar \tau^2
static double sim_kin_polybead_fw(struct ensem const* const ens,
    size_t const ipoly, size_t const ibead) {
  double const M = ens->R[ipoly].m / (2.0 * gsl_pow_2(ens->tau));

  if (ibead == ens->nmemb.bead - 1) {
    size_t const jpoly = ens->R[ipoly].ito;
    size_t const jbead = 0;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return M * sim_dist2(ens,
          &ens->R[ipoly].r[ibead], &ens->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead + 1;

    return M * sim_dist2(ens,
        &ens->R[ipoly].r[ibead], &ens->R[ipoly].r[jbead]);
  }
}

// Kinetic energy within a polymer from one bead to its two neighbors.
static double sim_kin_polybead(struct ensem const* const ens,
    size_t const ipoly, size_t const ibead) {
  return sim_kin_polybead_bw(ens, ipoly, ibead) +
    sim_kin_polybead_fw(ens, ipoly, ibead);
}

// Kinetic energy within a polymer.
static double sim_kin_poly(struct ensem const* const ens,
    size_t const ipoly) {
  double K = 0.0;

  for (size_t ibead = 0; ibead < ens->nmemb.bead; ++ibead)
    K += sim_kin_polybead_fw(ens, ipoly, ibead);

  return K;
}

// Kinetic energy.
static double sim_kin_total(struct ensem const* const ens) {
  double K = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    K += sim_kin_poly(ens, ipoly);

  return K;
}

// Potential energy between certain beads of all polymers.
static double sim_potint_bead(struct ensem const* const ens,
    size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < ens->nmemb.poly; ++jpoly)
      V += ens->Vint(ens,
          &ens->R[ipoly].r[ibead], &ens->R[jpoly].r[ibead]);

  return V;
}

// Additional potential energy between either of the ends of all polymers.
static double sim_potend_bead(struct ensem const* const ens,
    size_t const ibead) {
  double V = 0.0;

  if (ibead == 0) {
    for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
      if (ens->R[ipoly].ifrom == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < ens->nmemb.poly; ++jpoly)
          if (ens->R[jpoly].ifrom == SIZE_MAX)
            V += ens->Vend(ens,
                &ens->R[ipoly].r[ibead], &ens->R[jpoly].r[ibead]);
  } else if (ibead == ens->nmemb.bead - 1)
    for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
      if (ens->R[ipoly].ito == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < ens->nmemb.poly; ++jpoly)
          if (ens->R[jpoly].ito == SIZE_MAX)
            V += ens->Vend(ens,
                &ens->R[ipoly].r[ibead], &ens->R[jpoly].r[ibead]);

  return V;
}

// External potential energy for one bead in a polymer.
static double sim_potext_polybead(struct ensem const* const ens,
    size_t const ipoly, size_t const ibead) {
  return ens->Vext(ens, &ens->R[ipoly].r[ibead]);
}

// External potential energy for every bead in a polymer.
static double sim_potext_bead(struct ensem const* const ens, size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    V += sim_potext_polybead(ens, ipoly, ibead);

  return V;
}

// Potential energy for one bead over all polymers.
static double sim_pot_bead(struct ensem const* const ens,
    size_t const ibead) {
  return sim_potint_bead(ens, ibead) + sim_potend_bead(ens, ibead) +
    sim_potext_bead(ens, ibead);
}

// Potential energy.
static double sim_pot_total(struct ensem const* const ens) {
  double V = 0.0;

  for (size_t ibead = 0; ibead < ens->nmemb.bead; ++ibead)
    V += sim_pot_bead(ens, ibead);

  return V;
}

// Thermodynamic energy estimator for closed polymers.
// \langle E_T\rangle & = \frac 1{M \tau} \sum_{k = 1}^M
// \Bigl\langle\frac{d N} 2 - \frac{|R_{k + 1 \bmod M} - R_k|^2}{4 \lambda \tau} + \tau V(R_k)\Bigr\rangle
static double est_pimc_tde(struct ensem const* const ens) {
  double const E = (double) (ens->nmemb.dim * ens->nmemb.poly) /
    (2.0 * ens->tau);
  double const K = sim_kin_total(ens) /
    (double) (ens->nmemb.poly * ens->nmemb.bead);
  double const V = sim_pot_total(ens) /
    (double) (ens->nmemb.poly * ens->nmemb.bead);

  return E - K + V;
}

// TODO The crap thermodynamic energy estimator for open polymers.
static double est_pigs_crap(
    __attribute__ ((__unused__)) struct ensem const* const ens) {
  return NAN;
}

// TODO Thermodynamic energy estimator for open polymers.
static double est_pigs_tde(
    __attribute__ ((__unused__)) struct ensem const* const ens) {
  return NAN;
}

void sim_mass_const(struct sim* const sim, double const m) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
    sim->ens.R[ipoly].m = m;
}

void sim_perm_close(struct sim* const sim) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    sim->ens.R[ipoly].ifrom = ipoly;
    sim->ens.R[ipoly].ito = ipoly;
  }
}

void sim_perm_open(struct sim* const sim) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    sim->ens.R[ipoly].ifrom = SIZE_MAX;
    sim->ens.R[ipoly].ito = SIZE_MAX;
  }
}

void sim_perm_chain(struct sim* const sim) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    sim->ens.R[ipoly].ifrom = size_uwrap_dec(ipoly, sim->ens.nmemb.poly);
    sim->ens.R[ipoly].ito = size_uwrap_inc(ipoly, sim->ens.nmemb.poly);
  }
}

void sim_perm_random(struct sim* const sim) {
  size_t i[POLY_MAX];
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
    i[ipoly] = ipoly;

  gsl_ran_shuffle(sim->rng, i, sim->ens.nmemb.poly, sizeof *i);

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    sim->ens.R[i[ipoly]].ifrom = ipoly;
    sim->ens.R[ipoly].ito = i[ipoly];
  }
}

void sim_placer_point(struct sim* const sim, size_t const ipoly,
    struct bead const* const r, __attribute__ ((__unused__)) double const d,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      sim->ens.R[ipoly].r[ibead].d[idim] = r->d[idim];
}

void sim_placer_random(struct sim* const sim, size_t const ipoly,
    struct bead const* const r, double const d,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      sim->ens.R[ipoly].r[ibead].d[idim] = r->d[idim] +
        ran_open(sim->rng, d / 2.0);
}

void sim_placer_knot(struct sim* const sim, size_t const ipoly,
    struct bead const* const r, double const d, void const* const p) {
  double k[DIM_MAX];
  if (p == NULL)
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      k[idim] = 1.0;
  else
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      k[idim] = (double) ((size_t const*) p)[idim];

  for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead) {
    double const t = (double) ibead / (double) sim->ens.nmemb.bead;

    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim) {
      double const phi = (double) idim / (double) sim->ens.nmemb.dim;

      sim->ens.R[ipoly].r[ibead].d[idim] = r->d[idim] +
        d * cos(M_2PI * (k[idim] * t + phi)) / 2.0;
    }
  }
}

void sim_place_point(struct sim* const sim, sim_placer const f,
    void const* const p) {
  struct bead r;
  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
    r.d[idim] = 0.0;

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
    f(sim, ipoly, &r, sim->ens.L, p);
}

void sim_place_random(struct sim* const sim, sim_placer const f,
    void const* const p) {
  double const b = fp_rt(sim->ens.nmemb.poly, sim->ens.nmemb.dim);

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    struct bead r;
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      r.d[idim] = ran_open(sim->rng, sim->ens.L / 2.0);

    f(sim, ipoly, &r, sim->ens.L / b, p);
  }
}

void sim_place_lattice(struct sim* const sim, sim_placer const f,
    void const* const p) {
  size_t const nlin = size_cirt(sim->ens.nmemb.poly, sim->ens.nmemb.dim);

  double const a = -sim->ens.L / 2.0;
  double const b = (double) nlin;

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, nlin, i, sim->ens.nmemb.dim);

    struct bead r;
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim], 0.0, b, -a, a);

    f(sim, ipoly, &r, sim->ens.L / b, p);
  }
}

void sim_move_accept_ss(struct sim* const sim) {
  ++sim->params.ssm.acc;
}

void sim_move_reject_ss(struct sim* const sim) {
  size_t const ipoly = sim->hist.ssm.ipoly;
  size_t const ibead = sim->hist.ssm.ibead;

  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
    sim->ens.R[ipoly].r[ibead].d[idim] = sim->hist.ssm.r.d[idim];

  ++sim->params.ssm.rej;
}

void sim_move_adjust_ss(struct sim* const sim) {
  sim->params.ssm.h = fp_clamp(sim->params.ssm.h *
      fp_cbalance((double) sim->params.ssm.acc, (double) sim->params.ssm.rej),
      DBL_EPSILON, sim->ens.L);
}

void sim_move_ss(struct sim* const sim,
    size_t const ipoly, size_t const ibead) {
  sim->accept = sim_move_accept_ss;
  sim->reject = sim_move_reject_ss;
  sim->adjust = sim_move_adjust_ss;

  sim->hist.ssm.ipoly = ipoly;
  sim->hist.ssm.ibead = ibead;

  double (* const f)(double, double) =
    sim->ens.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim) {
    sim->hist.ssm.r.d[idim] = sim->ens.R[ipoly].r[ibead].d[idim];

    sim->ens.R[ipoly].r[ibead].d[idim] =
      f(sim->ens.R[ipoly].r[ibead].d[idim] +
          sim->params.ssm.h * ran_open(sim->rng, 1.0), sim->ens.L);
  }
}

void sim_move_accept_cmd(struct sim* const sim) {
  ++sim->params.cmd.acc;
}

void sim_move_reject_cmd(struct sim* const sim) {
  size_t const ipoly = sim->hist.cmd.ipoly;

  for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      sim->ens.R[ipoly].r[ibead].d[idim] =
        sim->hist.cmd.R.r[ibead].d[idim];

  ++sim->params.cmd.rej;
}

void sim_move_adjust_cmd(struct sim* const sim) {
  sim->params.cmd.h = fp_clamp(sim->params.cmd.h *
      fp_cbalance((double) sim->params.cmd.acc, (double) sim->params.cmd.rej),
      DBL_EPSILON, sim->ens.L);
}

void sim_move_cmd(struct sim* const sim, size_t const ipoly) {
  sim->accept = sim_move_accept_cmd;
  sim->reject = sim_move_reject_cmd;
  sim->adjust = sim_move_adjust_cmd;

  sim->hist.cmd.ipoly = ipoly;

  double (* const f)(double, double) =
    sim->ens.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim) {
    double const x = sim->params.cmd.h *
      ran_open(sim->rng, 1.0);

    for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead) {
      sim->hist.cmd.R.r[ibead].d[idim] =
        sim->ens.R[ipoly].r[ibead].d[idim];

      sim->ens.R[ipoly].r[ibead].d[idim] =
        f(sim->ens.R[ipoly].r[ibead].d[idim] + x, sim->ens.L);
    }
  }
}

void sim_move_accept_bisect(
    __attribute__ ((__unused__)) struct sim* const sim) {
  err_abort(NULL);
}

void sim_move_reject_bisect(
    __attribute__ ((__unused__)) struct sim* const sim) {
  err_abort(NULL);
}

void sim_move_bisect(
    __attribute__ ((__unused__)) struct sim* const sim,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

void sim_move_accept_swap(
    __attribute__ ((__unused__)) struct sim* const sim) {
  err_abort(NULL);
}

void sim_move_reject_swap(
    __attribute__ ((__unused__)) struct sim* const sim) {
  err_abort(NULL);
}

void sim_move_swap(
    __attribute__ ((__unused__)) struct sim* const sim,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

double sim_propose_ss(struct sim* const sim) {
  size_t const ipoly = ran_index(sim->rng, sim->ens.nmemb.poly);
  size_t const ibead = ran_index(sim->rng, sim->ens.nmemb.bead);

  double const V0 =
    sim_potint_bead(&sim->ens, ibead) + sim_potend_bead(&sim->ens, ibead) +
    sim_potext_polybead(&sim->ens, ipoly, ibead);
  double const K0 = sim_kin_polybead(&sim->ens, ipoly, ibead);

  sim_move_ss(sim, ipoly, ibead);

  double const V1 =
    sim_potint_bead(&sim->ens, ibead) + sim_potend_bead(&sim->ens, ibead) +
    sim_potext_polybead(&sim->ens, ipoly, ibead);
  double const K1 = sim_kin_polybead(&sim->ens, ipoly, ibead);

  return sim->ens.tau * (K1 + V1 - K0 - V0);
}

double sim_propose_cmd(struct sim* const sim) {
  size_t const ipoly = ran_index(sim->rng, sim->ens.nmemb.poly);

  double const V0 = sim_pot_total(&sim->ens);

  sim_move_cmd(sim, ipoly);

  double const V1 = sim_pot_total(&sim->ens);

  return sim->ens.tau * (V1 - V0);
}

static void decide(struct sim* const sim, double const DeltaS) {
  if (DeltaS > 0.0 && exp(-DeltaS) < ran_uopen(sim->rng, 1.0))
    sim->reject(sim);
  else
    sim->accept(sim);

  sim->adjust(sim);
}

static void posdist_accum(struct sim const* const sim) {
  double a;
  double b;
  ensem_extents(&sim->ens, &a, &b);

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead) {
      size_t i[DIM_MAX];

      for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
        i[idim] = size_uclamp(
            (size_t) (floor(fp_lerp(sim->ens.R[ipoly].r[ibead].d[idim],
                  a, b, 0.0, (double) sim->ens.nmemb.subdiv))),
            sim->ens.nmemb.subdiv);

      ++sim->p[size_unhc(sim->ens.nmemb.subdiv, i,
          sim->ens.nmemb.dim)];
    }
}

static void paircorr_accum(struct sim const* const sim) {
  size_t const nbin = size_pow(sim->ens.nmemb.subdiv, sim->ens.nmemb.dim);

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < sim->ens.nmemb.poly; ++jpoly)
      for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead) {
        size_t const i = size_uclamp(
            (size_t) (floor(fp_lerp(sim_dist(&sim->ens,
                    &sim->ens.R[ipoly].r[ibead],
                    &sim->ens.R[jpoly].r[ibead]),
                  0.0, (double) sim->ens.L, 0.0, (double) nbin))), nbin);

        ++sim->g[i];
      }
}

static void ensem_extents(struct ensem const* const ens,
    double* const a, double* const b) {
  if (ens->periodic) {
    *a = 0.0;
    *b = ens->L;
  } else {
    *a = -ens->L / 2.0;
    *b = ens->L / 2.0;
  }
}

bool res_close(__attribute__ ((__unused__)) struct sim const* const sim,
    FILE* const fp) {
  if (fclose(fp) == EOF)
    return false;

  return true;
}

FILE* res_open(struct sim const* const sim, char const* const str) {
  // Sanity check.
  if (strchr(str, '.') != NULL || strchr(str, '/') != NULL)
    return NULL;

  char buf[BUFSIZ];
  int const k = snprintf(buf, sizeof buf,
      "%s/%s.data", sim->id, str);
  if (k < 0 || (size_t) k >= sizeof buf)
    return NULL;

  return fopen(buf, "w");
}

bool res_disp(struct sim const* const sim, char const* const str,
    sim_printer const f, void const* const p) {
  FILE* const fp = res_open(sim, str);
  if (fp == NULL)
    return false;

  bool const q = f(sim, fp, p);

  if (!res_close(sim, fp))
    return false;

  return q;
}

bool print_periodic(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%d\n", sim->ens.periodic ? 1 : 0) < 0)
    return false;

  return true;
}

bool print_length(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%g\n", sim->ens.L) < 0)
    return false;

  return true;
}

bool print_pots(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  double a;
  double b;
  ensem_extents(&sim->ens, &a, &b);

  size_t const npt = size_pow(sim->ens.nmemb.subdiv, sim->ens.nmemb.dim);

  struct bead r0;
  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
    r0.d[idim] = 0.0;

  for (size_t ipt = 0; ipt < npt; ++ipt) {
    size_t i[DIM_MAX];
    size_hc(ipt, sim->ens.nmemb.subdiv, i, sim->ens.nmemb.dim);

    struct bead r;

    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim],
          0.0, (double) sim->ens.nmemb.subdiv, a, b);

    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      if (fprintf(fp, "%g ", r.d[idim]) < 0)
        return false;

    if (fprintf(fp, "%g %g\n",
          sim->ens.Vext(&sim->ens, &r),
          sim->ens.Vint(&sim->ens, &r0, &r)) < 0)
      return false;
  }

  return true;
}

bool print_ndim(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%zu\n", sim->ens.nmemb.dim) < 0)
    return false;

  return true;
}

bool print_npoly(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%zu\n", sim->ens.nmemb.poly) < 0)
    return false;

  return true;
}

bool print_nbead(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%zu\n", sim->ens.nmemb.bead) < 0)
    return false;

  return true;
}

bool print_nsubdiv(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%zu\n", sim->ens.nmemb.subdiv) < 0)
    return false;

  return true;
}

bool print_energy(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%zu %g %g %g\n",
        sim->ens.istep.prod, est_pimc_tde(&sim->ens),
        stats_mean(sim->tde), stats_sem(sim->tde)) < 0)
    return false;

  return true;
}

bool print_energy_corrtime(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%g\n", stats_corrtime(sim->tde)) < 0)
    return false;

  return true;
}

bool print_params(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%zu %zu %zu %zu %g %zu %zu %g\n",
        sim->ens.istep.thrm, sim->ens.istep.prod,
        sim->params.ssm.acc, sim->params.ssm.rej,
        sim->params.ssm.h,
        sim->params.cmd.acc, sim->params.cmd.rej,
        sim->params.cmd.h) < 0)
    return false;

  return true;
}

bool print_posdist(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  double a;
  double b;
  ensem_extents(&sim->ens, &a, &b);

  size_t const nbin = size_pow(sim->ens.nmemb.subdiv, sim->ens.nmemb.dim);

  size_t N = 0;
  for (size_t ibin = 0; ibin < nbin; ++ibin)
    N += sim->p[ibin];
  double const S = (N == 0 ? 1.0 : (double) N) * sim->ens.L;

  for (size_t ibin = 0; ibin < nbin; ++ibin) {
    size_t i[DIM_MAX];
    size_hc(ibin, sim->ens.nmemb.subdiv, i, sim->ens.nmemb.dim);

    struct bead r;

    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim] + 0.5,
          0.0, (double) sim->ens.nmemb.subdiv, a, b);

    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      if (fprintf(fp, "%g ", r.d[idim]) < 0)
        return false;

    if (fprintf(fp, "%g\n", (double) sim->p[ibin] / S) < 0)
      return false;
  }

  return true;
}

bool print_paircorr(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  size_t const nbin = size_pow(sim->ens.nmemb.subdiv, sim->ens.nmemb.dim);

  size_t N = 0.0;
  for (size_t ibin = 0; ibin < nbin; ++ibin)
    N += sim->g[ibin];
  double const S = N == 0 ? 1.0 : (double) N;

  for (size_t ibin = 0; ibin < nbin; ++ibin) {
    double const r = fp_lerp((double) ibin + 0.5,
        0.0, (double) nbin, 0.0, sim->ens.L);

    if (fprintf(fp, "%g %g\n", r, (double) sim->g[ibin] / S) < 0)
      return false;
  }

  return true;
}

static bool print_polys1(struct sim const* const sim, FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  if (fprintf(fp, "%g", (double) iindex * sim->ens.tau) < 0)
    return false;

  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
    if (fprintf(fp, " %g", sim->ens.R[ipoly].r[ibead].d[idim]) < 0)
      return false;

  if (fprintf(fp, "\n") < 0)
    return false;

  return true;
}

bool print_polys(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
      if (!print_polys1(sim, fp, ibead, ipoly, ibead))
        return false;

    if (sim->ens.R[ipoly].ito != SIZE_MAX)
      if (!print_polys1(sim, fp, sim->ens.nmemb.bead, sim->ens.R[ipoly].ito, 0))
        return false;

    if (fprintf(fp, "\n") < 0)
      return false;
  }

  return true;
}

bool print_progress(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  size_t const i = sim->ens.istep.thrm + sim->ens.istep.prod;
  size_t const n = sim->ens.nstep.thrm + sim->ens.nstep.prod;

  if (fprintf(fp, "(i_T + i_P) / (n_T + n_P) = "
        "(%zu + %zu) / (%zu + %zu) = %zu %%\n",
        sim->ens.istep.thrm, sim->ens.istep.prod,
        sim->ens.nstep.thrm, sim->ens.nstep.prod,
        n == 0 ? 100 : 100 * i / n) < 0)
    return false;

  return true;
}

bool print_results(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "E = %g +- %g (kappa = %g)\n",
        stats_mean(sim->tde), stats_corrsem(sim->tde),
        stats_corrtime(sim->tde)) < 0)
    return false;

  return true;
}

bool print_wrong_results_fast(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "E = %g +- %g (kappa = %g)\n",
        stats_mean(sim->tde), 100.0 * stats_sem(sim->tde), 100.0) < 0)
    return false;

  return true;
}

static bool prepare_save(struct sim const* const sim) {
  // Sanity check.
  if (strncmp(sim->id, "run", 3) != 0)
    err_abort(strncmp);

  if (symlink(sim->id, sim->id) == -1)
    err_abort(symlink);

  if (rename(sim->id, "run-latest") == -1)
    err_abort(rename);

  if (mkdir(sim->id, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    err_abort(mkdir);

  return true;
}

static bool save_const(struct sim const* const sim) {
  return res_disp(sim, "periodic", print_periodic, NULL) &&
    res_disp(sim, "length", print_length, NULL) &&
    res_disp(sim, "pots", print_pots, NULL) &&
    res_disp(sim, "ndim", print_ndim, NULL) &&
    res_disp(sim, "npoly", print_npoly, NULL) &&
    res_disp(sim, "nbead", print_nbead, NULL) &&
    res_disp(sim, "nsubdiv", print_nsubdiv, NULL);
}

static bool save_mut(struct sim const* const sim) {
  return res_disp(sim, "energy-corrtime", print_energy_corrtime, NULL) &&
    res_disp(sim, "posdist", print_posdist, NULL) &&
    res_disp(sim, "paircorr", print_paircorr, NULL) &&
    res_disp(sim, "polys", print_polys, NULL);
}

void sim_free(struct sim* const sim) {
  if (sim != NULL) {
    free(sim->g);
    free(sim->p);

    stats_free(sim->tde);

    gsl_rng_free(sim->rng);
  }

  free(sim);
}

struct sim* sim_alloc(size_t const ndim, size_t const npoly,
    size_t const nbead, size_t const nsubdiv,
    size_t const nthrm, size_t const nprod,
    size_t const nthrmrec, size_t const nprodrec,
    double const L, double const m, double const beta) {
  // These assertions guarantee that adding or multiplying two steps or indices
  // never wraps around, which makes their manipulation easier.
  dynamic_assert(ndim <= DIM_MAX, "too many dimensions");
  dynamic_assert(npoly <= POLY_MAX, "too many polymers");
  dynamic_assert(nbead <= BEAD_MAX, "too many beads");
  dynamic_assert(nthrm <= SQRT_SIZE_MAX, "too many thermalization steps");
  dynamic_assert(nprod <= SQRT_SIZE_MAX, "too many production steps");
  dynamic_assert(nthrmrec <= nthrm, "too many thermalization recording steps");
  dynamic_assert(nprodrec <= nprod, "too many production recording steps");
  dynamic_assert(L > 0.0, "empty simulation box");
  dynamic_assert(m > 0.0, "negative default mass");
  dynamic_assert(beta > 0.0, "weird things going on");

  bool p = true;

  struct sim* const sim = malloc(sizeof *sim);
  if (sim == NULL)
    p = false;
  else {
    sim->rng = gsl_rng_alloc(gsl_rng_env_setup());
    if (sim->rng == NULL)
      p = false;
    else if (!ran_dateid(sim->rng, "run", sim->id, sizeof sim->id))
      p = false;

    sim->accept = NULL;
    sim->reject = NULL;
    sim->adjust = NULL;

    sim->tde = stats_alloc(nprod);
    if (sim->tde == NULL)
      p = false;

    size_t const nbin = size_pow(nsubdiv, ndim);
    sim->p = malloc(nbin * sizeof *sim->p);
    if (sim->p == NULL)
      p = false;
    else
      for (size_t ibin = 0; ibin < nbin; ++ibin)
        sim->p[ibin] = 0.0;

    sim->g = malloc(nbin * sizeof *sim->g);
    if (sim->g == NULL)
      p = false;
    else
      for (size_t ibin = 0; ibin < nbin; ++ibin)
        sim->g[ibin] = 0.0;

    sim->params.ssm.acc = 0;
    sim->params.ssm.rej = 0;
    sim->params.ssm.h = L / 2.0;

    sim->params.cmd.acc = 0;
    sim->params.cmd.rej = 0;
    sim->params.cmd.h = L / 2.0;

    sim->ens.periodic = false;
    sim->ens.L = L;
    sim->ens.tau = beta / (double) nbead;

    sim->ens.Vint = sim_pot_zero;
    sim->ens.Vend = sim_pot_zero;
    sim->ens.Vext = sim_potext_zero;

    sim->ens.nmemb.dim = ndim;
    sim->ens.nmemb.poly = npoly;
    sim->ens.nmemb.bead = nbead;
    sim->ens.nmemb.subdiv = nsubdiv;

    sim->ens.nstep.thrm = nthrm;
    sim->ens.nstep.prod = nprod;
    sim->ens.nstep.thrmrec = nthrmrec;
    sim->ens.nstep.prodrec = nprodrec;

    sim->ens.istep.thrm = 0;
    sim->ens.istep.prod = 0;
    sim->ens.istep.thrmrec = 0;
    sim->ens.istep.prodrec = 0;

    sim_mass_const(sim, m);
    sim_perm_close(sim);
    sim_place_point(sim, sim_placer_point, NULL);
  }

  if (p)
    return sim;
  else {
    sim_free(sim);

    return NULL;
  }
}

bool sim_run(struct sim* const sim) {
  err_reset();

#ifdef _GNU_SOURCE
#ifdef DEBUG
  int const excepts = feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if (excepts == -1)
    err_abort(feenableexcept);
#endif
#endif

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

  if (!prepare_save(sim))
    err_abort(prepare_save);

  if (!save_const(sim))
    err_abort(save_const);

  FILE* const paramsfp = res_open(sim, "params");
  if (paramsfp == NULL)
    err_abort(res_open);

  FILE* const energyfp = res_open(sim, "energy");
  if (energyfp == NULL)
    err_abort(res_open);

  sim_proposer const proposers[] = {sim_propose_ss, sim_propose_cmd};
  size_t const nproposer = sizeof proposers / sizeof *proposers;

  for (size_t istep = 0;
      istep < sim->ens.nstep.thrm + sim->ens.nstep.prod;
      ++istep) {
    decide(sim, proposers[ran_index(sim->rng, nproposer)](sim));

    int signum;
    if (sigs_use(&signum))
      switch (signum) {
        case SIGUSR1:
          (void) print_progress(sim, stdout, NULL);
          break;
        case SIGUSR2:
          (void) save_mut(sim);
          break;
      }

    if (istep < sim->ens.nstep.thrm) {
      if (sim->ens.nstep.prodrec * sim->ens.istep.thrm >
          sim->ens.nstep.prod * sim->ens.istep.thrmrec) {
        if (!print_params(sim, paramsfp, NULL))
          err_abort(print_params);

        ++sim->ens.istep.thrmrec;
      }

      ++sim->ens.istep.thrm;
    } else {
      (void) stats_accum(sim->tde, est_pimc_tde(&sim->ens));

      posdist_accum(sim);
      paircorr_accum(sim);

      if (sim->ens.nstep.prodrec * sim->ens.istep.prod >
          sim->ens.nstep.prod * sim->ens.istep.prodrec) {
        if (!print_energy(sim, energyfp, NULL))
          err_abort(print_energy);

        if (!print_params(sim, paramsfp, NULL))
          err_abort(print_params);

        ++sim->ens.istep.prodrec;
      }

      ++sim->ens.istep.prod;
    }
  }

  if (!res_close(sim, energyfp))
    err_abort(res_close);

  if (!res_close(sim, paramsfp))
    err_abort(res_close);

  if (!save_mut(sim))
    err_abort(save_mut);

  if (!print_progress(sim, stdout, NULL))
    err_abort(print_progress);

  if (!print_wrong_results_fast(sim, stdout, NULL))
    err_abort(print_wrong_results_fast);

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    err_abort(feenableexcept);
#endif
#endif

  return true;
}
