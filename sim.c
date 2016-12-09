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

// This structure represents an ensemble and is a part of `struct napkin`.
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
struct napkin {
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

double bead_norm2(struct ensem const* const ens,
    struct bead const* const r) {
  double s = 0.0;

  double (* const f)(double, double) =
    ens->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ens->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r->d[idim], ens->L));

  return s;
}

double bead_norm(struct ensem const* const ens,
    struct bead const* const r) {
  return sqrt(bead_norm2(ens, r));
}

double bead_dist2(struct ensem const* const ens,
    struct bead const* const r0, struct bead const* const r1) {
  double s = 0.0;

  double (* const f)(double, double) =
    ens->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ens->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r1->d[idim] - r0->d[idim], ens->L));

  return s;
}

double bead_dist(struct ensem const* const ens,
    struct bead const* const r0, struct bead const* const r1) {
  return sqrt(bead_dist2(ens, r0, r1));
}

struct ensem* sim_get_ensem(struct napkin* const nap) {
  return &nap->ens;
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
      return M * bead_dist2(ens,
          &ens->R[ipoly].r[ibead], &ens->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead - 1;

    return M * bead_dist2(ens,
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
      return M * bead_dist2(ens,
          &ens->R[ipoly].r[ibead], &ens->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead + 1;

    return M * bead_dist2(ens,
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
static double est_pigs_crap(struct ensem const* const ens) {
  return NAN;
}

// TODO Thermodynamic energy estimator for open polymers.
static double est_pigs_tde(struct ensem const* const ens) {
  return NAN;
}

void sim_mass_const(struct napkin* const nap, double const m) {
  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly)
    nap->ens.R[ipoly].m = m;
}

void sim_perm_close(struct napkin* const nap) {
  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly) {
    nap->ens.R[ipoly].ifrom = ipoly;
    nap->ens.R[ipoly].ito = ipoly;
  }
}

void sim_perm_open(struct napkin* const nap) {
  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly) {
    nap->ens.R[ipoly].ifrom = SIZE_MAX;
    nap->ens.R[ipoly].ito = SIZE_MAX;
  }
}

void sim_perm_chain(struct napkin* const nap) {
  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly) {
    nap->ens.R[ipoly].ifrom = size_uwrap_dec(ipoly, nap->ens.nmemb.poly);
    nap->ens.R[ipoly].ito = size_uwrap_inc(ipoly, nap->ens.nmemb.poly);
  }
}

void sim_perm_random(struct napkin* const nap) {
  size_t i[POLY_MAX];
  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly)
    i[ipoly] = ipoly;

  gsl_ran_shuffle(nap->rng, i, nap->ens.nmemb.poly, sizeof *i);

  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly) {
    nap->ens.R[i[ipoly]].ifrom = ipoly;
    nap->ens.R[ipoly].ito = i[ipoly];
  }
}

void sim_placer_point(struct napkin* const nap, size_t const ipoly,
    struct bead const* const r, __attribute__ ((__unused__)) double const d) {
  for (size_t ibead = 0; ibead < nap->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      nap->ens.R[ipoly].r[ibead].d[idim] = r->d[idim];
}

void sim_placer_random(struct napkin* const nap, size_t const ipoly,
    struct bead const* const r, double const d) {
  for (size_t ibead = 0; ibead < nap->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      nap->ens.R[ipoly].r[ibead].d[idim] = r->d[idim] +
        ran_open(nap->rng, d / 2.0);
}

void sim_placer_unknot(struct napkin* const nap, size_t const ipoly,
    struct bead const* const r, double const d) {
  for (size_t ibead = 0; ibead < nap->ens.nmemb.bead; ++ibead) {
    double const t = (double) ibead / (double) nap->ens.nmemb.bead;

    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim) {
      double const phi = (double) idim / (double) nap->ens.nmemb.dim;

      nap->ens.R[ipoly].r[ibead].d[idim] = r->d[idim] +
        (d / 2.0) * sin(M_2PI * (t + phi / 2.0));
    }
  }
}

void sim_placer_tinyuk(struct napkin* const nap, size_t const ipoly,
    struct bead const* const r, double const d) {
  sim_placer_unknot(nap, ipoly, r, d / 2.0);
}

void sim_place_point(struct napkin* const nap,
    void (* const f)(struct napkin*, size_t, struct bead const*, double)) {
  struct bead r;
  for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
    r.d[idim] = 0.0;

  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly)
    f(nap, ipoly, &r, nap->ens.L);
}

void sim_place_random(struct napkin* const nap,
    void (* const f)(struct napkin*, size_t, struct bead const*, double)) {
  double const b = fp_rt(nap->ens.nmemb.poly, nap->ens.nmemb.dim);

  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly) {
    struct bead r;
    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      r.d[idim] = ran_open(nap->rng, nap->ens.L / 2.0);

    f(nap, ipoly, &r, nap->ens.L / b);
  }
}

void sim_place_lattice(struct napkin* const nap,
    void (* const f)(struct napkin*, size_t, struct bead const*, double)) {
  size_t const nlin = size_cirt(nap->ens.nmemb.poly, nap->ens.nmemb.dim);

  double const a = -nap->ens.L / 2.0;
  double const b = (double) nlin;

  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, nlin, i, nap->ens.nmemb.dim);

    struct bead r;
    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim], 0.0, b, -a, a);

    f(nap, ipoly, &r, nap->ens.L / b);
  }
}

void sim_move_accept_ss(struct napkin* const nap) {
  ++nap->params.ssm.acc;
}

void sim_move_reject_ss(struct napkin* const nap) {
  size_t const ipoly = nap->hist.ssm.ipoly;
  size_t const ibead = nap->hist.ssm.ibead;

  for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
    nap->ens.R[ipoly].r[ibead].d[idim] = nap->hist.ssm.r.d[idim];

  ++nap->params.ssm.rej;
}

void sim_move_adjust_ss(struct napkin* const nap) {
  nap->params.ssm.h = fp_clamp(nap->params.ssm.h *
      fp_cbalance((double) nap->params.ssm.acc, (double) nap->params.ssm.rej),
      DBL_EPSILON, nap->ens.L);
}

void sim_move_ss(struct napkin* const nap,
    size_t const ipoly, size_t const ibead) {
  nap->accept = sim_move_accept_ss;
  nap->reject = sim_move_reject_ss;
  nap->adjust = sim_move_adjust_ss;

  nap->hist.ssm.ipoly = ipoly;
  nap->hist.ssm.ibead = ibead;

  double (* const f)(double, double) =
    nap->ens.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim) {
    nap->hist.ssm.r.d[idim] = nap->ens.R[ipoly].r[ibead].d[idim];

    nap->ens.R[ipoly].r[ibead].d[idim] =
      f(nap->ens.R[ipoly].r[ibead].d[idim] +
          nap->params.ssm.h * ran_open(nap->rng, 1.0), nap->ens.L);
  }
}

void sim_move_accept_cmd(struct napkin* const nap) {
  ++nap->params.cmd.acc;
}

void sim_move_reject_cmd(struct napkin* const nap) {
  size_t const ipoly = nap->hist.cmd.ipoly;

  for (size_t ibead = 0; ibead < nap->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      nap->ens.R[ipoly].r[ibead].d[idim] =
        nap->hist.cmd.R.r[ibead].d[idim];

  ++nap->params.cmd.rej;
}

void sim_move_adjust_cmd(struct napkin* const nap) {
  nap->params.cmd.h = fp_clamp(nap->params.cmd.h *
      fp_cbalance((double) nap->params.cmd.acc, (double) nap->params.cmd.rej),
      DBL_EPSILON, nap->ens.L);
}

void sim_move_cmd(struct napkin* const nap, size_t const ipoly) {
  nap->accept = sim_move_accept_cmd;
  nap->reject = sim_move_reject_cmd;
  nap->adjust = sim_move_adjust_cmd;

  nap->hist.cmd.ipoly = ipoly;

  double (* const f)(double, double) =
    nap->ens.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim) {
    double const x = nap->params.cmd.h *
      ran_open(nap->rng, 1.0);

    for (size_t ibead = 0; ibead < nap->ens.nmemb.bead; ++ibead) {
      nap->hist.cmd.R.r[ibead].d[idim] =
        nap->ens.R[ipoly].r[ibead].d[idim];

      nap->ens.R[ipoly].r[ibead].d[idim] =
        f(nap->ens.R[ipoly].r[ibead].d[idim] + x, nap->ens.L);
    }
  }
}

void sim_move_accept_bisect(
    __attribute__ ((__unused__)) struct napkin* const nap) {
  err_abort(NULL);
}

void sim_move_reject_bisect(
    __attribute__ ((__unused__)) struct napkin* const nap) {
  err_abort(NULL);
}

void sim_move_bisect(
    __attribute__ ((__unused__)) struct napkin* const nap,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

void sim_move_accept_swap(
    __attribute__ ((__unused__)) struct napkin* const nap) {
  err_abort(NULL);
}

void sim_move_reject_swap(
    __attribute__ ((__unused__)) struct napkin* const nap) {
  err_abort(NULL);
}

void sim_move_swap(
    __attribute__ ((__unused__)) struct napkin* const nap,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

double sim_propose_ss(struct napkin* const nap) {
  size_t const ipoly = ran_index(nap->rng, nap->ens.nmemb.poly);
  size_t const ibead = ran_index(nap->rng, nap->ens.nmemb.bead);

  double const V0 =
    sim_potint_bead(&nap->ens, ibead) + sim_potend_bead(&nap->ens, ibead) +
    sim_potext_polybead(&nap->ens, ipoly, ibead);
  double const K0 = sim_kin_polybead(&nap->ens, ipoly, ibead);

  sim_move_ss(nap, ipoly, ibead);

  double const V1 =
    sim_potint_bead(&nap->ens, ibead) + sim_potend_bead(&nap->ens, ibead) +
    sim_potext_polybead(&nap->ens, ipoly, ibead);
  double const K1 = sim_kin_polybead(&nap->ens, ipoly, ibead);

  return nap->ens.tau * (K1 + V1 - K0 - V0);
}

double sim_propose_cmd(struct napkin* const nap) {
  size_t const ipoly = ran_index(nap->rng, nap->ens.nmemb.poly);

  double const V0 = sim_pot_total(&nap->ens);

  sim_move_cmd(nap, ipoly);

  double const V1 = sim_pot_total(&nap->ens);

  return nap->ens.tau * (V1 - V0);
}

static void decide(struct napkin* const nap, double const DeltaS) {
  if (DeltaS > 0.0 && exp(-DeltaS) < ran_uopen(nap->rng, 1.0))
    nap->reject(nap);
  else
    nap->accept(nap);

  nap->adjust(nap);
}

static void posdist_accum(struct napkin const* const nap) {
  double a;
  double b;
  ensem_extents(&nap->ens, &a, &b);

  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < nap->ens.nmemb.bead; ++ibead) {
      size_t i[DIM_MAX];

      for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
        i[idim] = size_uclamp(
            (size_t) (floor(fp_lerp(nap->ens.R[ipoly].r[ibead].d[idim],
                  a, b, 0.0, (double) nap->ens.nmemb.subdiv))),
            nap->ens.nmemb.subdiv);

      ++nap->p[size_unhc(nap->ens.nmemb.subdiv, i,
          nap->ens.nmemb.dim)];
    }
}

static void paircorr_accum(struct napkin const* const nap) {
  size_t const nbin = size_pow(nap->ens.nmemb.subdiv,
      nap->ens.nmemb.dim);

  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < nap->ens.nmemb.poly; ++jpoly)
      for (size_t ibead = 0; ibead < nap->ens.nmemb.bead; ++ibead) {
        size_t const i = size_uclamp(
            (size_t) (floor(fp_lerp(bead_dist(&nap->ens,
                    &nap->ens.R[ipoly].r[ibead],
                    &nap->ens.R[jpoly].r[ibead]),
                  0.0, (double) nap->ens.L, 0.0, (double) nbin))), nbin);

        ++nap->g[i];
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

static bool res_close(
    __attribute__ ((__unused__)) struct napkin const* const nap,
    FILE* const fp) {
  if (fclose(fp) == EOF)
    return false;

  return true;
}

static FILE* res_open(struct napkin const* const nap,
    char const* const str) {
  // Sanity check.
  if (strchr(str, '.') != NULL || strchr(str, '/') != NULL)
    return NULL;

  char buf[BUFSIZ];
  int const k = snprintf(buf, sizeof buf,
      "%s/%s.data", nap->id, str);
  if (k < 0 || (size_t) k >= sizeof buf)
    return NULL;

  return fopen(buf, "w");
}

static bool res_disp(struct napkin const* const nap,
    char const* const str, sim_printer const f) {
  FILE* const fp = res_open(nap, str);
  if (fp == NULL)
    return false;

  bool const p = f(nap, fp);

  if (!res_close(nap, fp))
    return false;

  return p;
}

static bool print_periodic(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "%d\n", nap->ens.periodic ? 1 : 0) < 0)
    return false;

  return true;
}

static bool print_length(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "%g\n", nap->ens.L) < 0)
    return false;

  return true;
}

static bool print_pots(struct napkin const* const nap, FILE* const fp) {
  double a;
  double b;
  ensem_extents(&nap->ens, &a, &b);

  size_t const npt = size_pow(nap->ens.nmemb.subdiv, nap->ens.nmemb.dim);

  struct bead r0;
  for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
    r0.d[idim] = 0.0;

  for (size_t ipt = 0; ipt < npt; ++ipt) {
    size_t i[DIM_MAX];
    size_hc(ipt, nap->ens.nmemb.subdiv, i, nap->ens.nmemb.dim);

    struct bead r;

    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim],
          0.0, (double) nap->ens.nmemb.subdiv, a, b);

    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      if (fprintf(fp, "%g ", r.d[idim]) < 0)
        return false;

    if (fprintf(fp, "%g %g\n",
          nap->ens.Vext(&nap->ens, &r),
          nap->ens.Vint(&nap->ens, &r0, &r)) < 0)
      return false;
  }

  return true;
}

static bool print_ndim(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "%zu\n", nap->ens.nmemb.dim) < 0)
    return false;

  return true;
}

static bool print_npoly(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "%zu\n", nap->ens.nmemb.poly) < 0)
    return false;

  return true;
}

static bool print_nbead(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "%zu\n", nap->ens.nmemb.bead) < 0)
    return false;

  return true;
}

static bool print_nsubdiv(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "%zu\n", nap->ens.nmemb.subdiv) < 0)
    return false;

  return true;
}

static bool print_energy(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "%zu %g %g %g\n",
        nap->ens.istep.prod, est_pimc_tde(&nap->ens),
        stats_mean(nap->tde), stats_sem(nap->tde)) < 0)
    return false;

  return true;
}

static bool print_energy_corrtime(struct napkin const* const nap,
    FILE* const fp) {
  if (fprintf(fp, "%g\n", stats_corrtime(nap->tde)) < 0)
    return false;

  return true;
}

static bool print_params(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "%zu %zu %zu %zu %g %zu %zu %g\n",
        nap->ens.istep.thrm, nap->ens.istep.prod,
        nap->params.ssm.acc, nap->params.ssm.rej,
        nap->params.ssm.h,
        nap->params.cmd.acc, nap->params.cmd.rej,
        nap->params.cmd.h) < 0)
    return false;

  return true;
}

static bool print_posdist(struct napkin const* const nap, FILE* const fp) {
  double a;
  double b;
  ensem_extents(&nap->ens, &a, &b);

  size_t const nbin = size_pow(nap->ens.nmemb.subdiv, nap->ens.nmemb.dim);

  size_t N = 0;
  for (size_t ibin = 0; ibin < nbin; ++ibin)
    N += nap->p[ibin];
  double const S = (N == 0 ? 1.0 : (double) N) * nap->ens.L;

  for (size_t ibin = 0; ibin < nbin; ++ibin) {
    size_t i[DIM_MAX];
    size_hc(ibin, nap->ens.nmemb.subdiv, i, nap->ens.nmemb.dim);

    struct bead r;

    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim] + 0.5,
          0.0, (double) nap->ens.nmemb.subdiv, a, b);

    for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
      if (fprintf(fp, "%g ", r.d[idim]) < 0)
        return false;

    if (fprintf(fp, "%g\n", (double) nap->p[ibin] / S) < 0)
      return false;
  }

  return true;
}

static bool print_paircorr(struct napkin const* const nap, FILE* const fp) {
  size_t const nbin = size_pow(nap->ens.nmemb.subdiv, nap->ens.nmemb.dim);

  size_t N = 0.0;
  for (size_t ibin = 0; ibin < nbin; ++ibin)
    N += nap->g[ibin];
  double const S = N == 0 ? 1.0 : (double) N;

  for (size_t ibin = 0; ibin < nbin; ++ibin) {
    double const r = fp_lerp((double) ibin + 0.5,
        0.0, (double) nbin, 0.0, nap->ens.L);

    if (fprintf(fp, "%g %g\n", r, (double) nap->g[ibin] / S) < 0)
      return false;
  }

  return true;
}

static bool print_polys1(struct napkin const* const nap, FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  if (fprintf(fp, "%g", (double) iindex * nap->ens.tau) < 0)
    return false;

  for (size_t idim = 0; idim < nap->ens.nmemb.dim; ++idim)
    if (fprintf(fp, " %g", nap->ens.R[ipoly].r[ibead].d[idim]) < 0)
      return false;

  if (fprintf(fp, "\n") < 0)
    return false;

  return true;
}

static bool print_polys(struct napkin const* const nap, FILE* const fp) {
  for (size_t ipoly = 0; ipoly < nap->ens.nmemb.poly; ++ipoly) {
    for (size_t ibead = 0; ibead < nap->ens.nmemb.bead; ++ibead)
      if (!print_polys1(nap, fp, ibead, ipoly, ibead))
        return false;

    if (nap->ens.R[ipoly].ito != SIZE_MAX)
      if (!print_polys1(nap, fp, nap->ens.nmemb.bead, nap->ens.R[ipoly].ito, 0))
        return false;

    if (fprintf(fp, "\n") < 0)
      return false;
  }

  return true;
}

static bool print_progress(struct napkin const* const nap, FILE* const fp) {
  size_t const i = nap->ens.istep.thrm + nap->ens.istep.prod;
  size_t const n = nap->ens.nstep.thrm + nap->ens.nstep.prod;

  if (fprintf(fp, "i / n = (%zu + %zu) / (%zu + %zu) = %zu %%\n",
        nap->ens.istep.thrm, nap->ens.istep.prod,
        nap->ens.nstep.thrm, nap->ens.nstep.prod,
        n == 0 ? 100 : 100 * i / n) < 0)
    return false;

  return true;
}

static bool print_results(struct napkin const* const nap, FILE* const fp) {
  if (fprintf(fp, "E = %g +- %g (kappa = %g)\n",
        stats_mean(nap->tde), stats_corrsem(nap->tde),
        stats_corrtime(nap->tde)) < 0)
    return false;

  return true;
}

static bool print_wrong_results_fast(struct napkin const* const nap,
    FILE* const fp) {
  if (fprintf(fp, "E = %g +- %g (kappa = %g)\n",
        stats_mean(nap->tde), 100.0 * stats_sem(nap->tde), 100.0) < 0)
    return false;

  return true;
}

static bool prepare_save(struct napkin const* const nap) {
  // Sanity check.
  if (strncmp(nap->id, "run", 3) != 0)
    err_abort(strncmp);

  if (symlink(nap->id, nap->id) == -1)
    err_abort(symlink);

  if (rename(nap->id, "run-latest") == -1)
    err_abort(rename);

  if (mkdir(nap->id, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    err_abort(mkdir);

  return true;
}

static bool save_const(struct napkin const* const nap) {
  return res_disp(nap, "periodic", print_periodic) &&
    res_disp(nap, "length", print_length) &&
    res_disp(nap, "pots", print_pots) &&
    res_disp(nap, "ndim", print_ndim) &&
    res_disp(nap, "npoly", print_npoly) &&
    res_disp(nap, "nbead", print_nbead) &&
    res_disp(nap, "nsubdiv", print_nsubdiv);
}

static bool save_mut(struct napkin const* const nap) {
  return res_disp(nap, "energy-corrtime", print_energy_corrtime) &&
    res_disp(nap, "posdist", print_posdist) &&
    res_disp(nap, "paircorr", print_paircorr) &&
    res_disp(nap, "polys", print_polys);
}

void napkin_free(struct napkin* const nap) {
  if (nap != NULL) {
    free(nap->g);
    free(nap->p);

    stats_free(nap->tde);

    gsl_rng_free(nap->rng);
  }

  free(nap);
}

struct napkin* napkin_alloc(size_t const ndim, size_t const npoly,
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

  struct napkin* const nap = malloc(sizeof *nap);
  if (nap == NULL)
    p = false;
  else {
    nap->rng = gsl_rng_alloc(gsl_rng_env_setup());
    if (nap->rng == NULL)
      p = false;
    else if (!ran_dateid(nap->rng, "run", nap->id, sizeof nap->id))
      p = false;

    nap->accept = NULL;
    nap->reject = NULL;
    nap->adjust = NULL;

    nap->tde = stats_alloc(nprod);
    if (nap->tde == NULL)
      p = false;

    size_t const nbin = size_pow(nsubdiv, ndim);
    nap->p = malloc(nbin * sizeof *nap->p);
    if (nap->p == NULL)
      p = false;
    else
      for (size_t ibin = 0; ibin < nbin; ++ibin)
        nap->p[ibin] = 0.0;

    nap->g = malloc(nbin * sizeof *nap->g);
    if (nap->g == NULL)
      p = false;
    else
      for (size_t ibin = 0; ibin < nbin; ++ibin)
        nap->g[ibin] = 0.0;

    nap->params.ssm.acc = 0;
    nap->params.ssm.rej = 0;
    nap->params.ssm.h = L / 2.0;

    nap->params.cmd.acc = 0;
    nap->params.cmd.rej = 0;
    nap->params.cmd.h = L / 2.0;

    nap->ens.periodic = false;
    nap->ens.L = L;
    nap->ens.tau = beta / (double) nbead;

    nap->ens.Vint = sim_pot_zero;
    nap->ens.Vend = sim_pot_zero;
    nap->ens.Vext = sim_potext_zero;

    nap->ens.nmemb.dim = ndim;
    nap->ens.nmemb.poly = npoly;
    nap->ens.nmemb.bead = nbead;
    nap->ens.nmemb.subdiv = nsubdiv;

    nap->ens.nstep.thrm = nthrm;
    nap->ens.nstep.prod = nprod;
    nap->ens.nstep.thrmrec = nthrmrec;
    nap->ens.nstep.prodrec = nprodrec;

    nap->ens.istep.thrm = 0;
    nap->ens.istep.prod = 0;
    nap->ens.istep.thrmrec = 0;
    nap->ens.istep.prodrec = 0;

    sim_mass_const(nap, m);
    sim_perm_close(nap);
    sim_place_point(nap, sim_placer_point);
  }

  if (p)
    return nap;
  else {
    napkin_free(nap);

    return NULL;
  }
}

bool sim_run(struct napkin* const nap) {
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

  /*
  if (!prepare_save(nap))
    err_abort(prepare_save);

  if (!save_const(nap))
    err_abort(save_const);
  */
  // TODO Delete this.
  strcpy(&nap->id, "run-2016");
  (void) mkdir("run-2016", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  nap->ens.L = 10.0;
  nap->ens.periodic = 1;
  sim_place_point(nap, sim_placer_point);
  sim_place_point(nap, sim_placer_unknot);
  sim_place_point(nap, sim_placer_random);
  sim_place_point(nap, sim_placer_tinyuk);
  sim_place_lattice(nap, sim_placer_point);
  sim_place_lattice(nap, sim_placer_random);
  sim_place_lattice(nap, sim_placer_unknot);
  sim_place_lattice(nap, sim_placer_tinyuk);
  sim_place_random(nap, sim_placer_point);
  sim_place_random(nap, sim_placer_random);
  sim_place_random(nap, sim_placer_unknot);
  sim_place_random(nap, sim_placer_tinyuk);
  sim_perm_random(nap);
  (void) res_disp(nap, "polys", print_polys);
  (void) save_const(nap);
  exit(0);

  FILE* const paramsfp = res_open(nap, "params");
  if (paramsfp == NULL)
    err_abort(res_open);

  FILE* const energyfp = res_open(nap, "energy");
  if (energyfp == NULL)
    err_abort(res_open);

  sim_proposer const proposers[] = {sim_propose_ss, sim_propose_cmd};
  size_t const nproposer = sizeof proposers / sizeof *proposers;

  for (size_t istep = 0;
      istep < nap->ens.nstep.thrm + nap->ens.nstep.prod;
      ++istep) {
    decide(nap, proposers[ran_index(nap->rng, nproposer)](nap));

    int signum;
    if (sigs_use(&signum))
      switch (signum) {
        case SIGUSR1:
          (void) print_progress(nap, stdout);
          break;
        case SIGUSR2:
          (void) save_mut(nap);
          break;
      }

    if (istep < nap->ens.nstep.thrm) {
      if (nap->ens.nstep.prodrec * nap->ens.istep.thrm >
          nap->ens.nstep.prod * nap->ens.istep.thrmrec) {
        if (!print_params(nap, paramsfp))
          err_abort(print_params);

        ++nap->ens.istep.thrmrec;
      }

      ++nap->ens.istep.thrm;
    } else {
      (void) stats_accum(nap->tde, est_pimc_tde(&nap->ens));

      posdist_accum(nap);
      paircorr_accum(nap);

      if (nap->ens.nstep.prodrec * nap->ens.istep.prod >
          nap->ens.nstep.prod * nap->ens.istep.prodrec) {
        if (!print_energy(nap, energyfp))
          err_abort(print_energy);

        if (!print_params(nap, paramsfp))
          err_abort(print_params);

        ++nap->ens.istep.prodrec;
      }

      ++nap->ens.istep.prod;
    }
  }

  if (!res_close(nap, energyfp))
    err_abort(res_close);

  if (!res_close(nap, paramsfp))
    err_abort(res_close);

  if (!save_mut(nap))
    err_abort(save_mut);

  if (!print_progress(nap, stdout))
    err_abort(print_progress);

  if (!print_wrong_results_fast(nap, stdout))
    err_abort(print_wrong_results_fast);

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    err_abort(feenableexcept);
#endif
#endif

  return true;
}
