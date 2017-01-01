#include "err.h"
#include "exts.h"
#include "fp.h"
#include "hist.h"
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
// `DIM_MAX`, `POLY_MAX`, `BEAD_MAX`, `SUBDIV_MAX` and `DIV_MAX`.
struct memb {
  size_t dim;
  size_t poly;
  size_t bead;
  size_t subdiv;
  size_t div;
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

struct bead {
  double r[DIM_MAX];
};

struct poly {
  double mass;
  size_t ifrom;
  size_t ito;
  struct bead r[BEAD_MAX];
};

struct ens {
  bool periodic;
  double length;
  double a;
  double b;
  double epsilon;
  ens_pot potint;
  ens_potext potext;
  struct memb nmemb;
  struct step nstep;
  struct step istep;
  struct poly r[POLY_MAX];
};

struct sim {
  gsl_rng* rng;
  struct ens ens;
  struct stats* tde;
  struct stats* tdei[BEAD_MAX];
  struct stats* mixed[BEAD_MAX];
  struct hist* p;
  struct hist* g;
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
      struct poly r;
    } cmd;
  } cache;
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
  char id[BUFSIZ];
};

typedef double (* sim_wrapper)(double, double);

double ens_pot_zero(
    __attribute__ ((__unused__)) struct ens const* const ens,
    __attribute__ ((__unused__)) struct bead const* const r0,
    __attribute__ ((__unused__)) struct bead const* const r1) {
  return 0.0;
}

double ens_potext_zero(
    __attribute__ ((__unused__)) struct ens const* const ens,
    __attribute__ ((__unused__)) struct bead const* const r) {
  return 0.0;
}

double ens_norm2(struct ens const* const ens, struct bead const* const r) {
  double s = 0.0;

  sim_wrapper const f = ens->periodic ? fp_swrap : fp_constant;

  for (size_t idim = 0; idim < ens->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r->r[idim], ens->length));

  return s;
}

double ens_norm(struct ens const* const ens, struct bead const* const r) {
  return sqrt(ens_norm2(ens, r));
}

double ens_dist2(struct ens const* const ens,
    struct bead const* const r0, struct bead const* const r1) {
  double s = 0.0;

  sim_wrapper const f = ens->periodic ? fp_swrap : fp_constant;

  for (size_t idim = 0; idim < ens->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r1->r[idim] - r0->r[idim], ens->length));

  return s;
}

double ens_dist(struct ens const* const ens,
    struct bead const* const r0, struct bead const* const r1) {
  return sqrt(ens_dist2(ens, r0, r1));
}

void sim_set_potint(struct sim* const sim, ens_pot const v) {
  sim->ens.potint = v;
}

void sim_set_potext(struct sim* const sim, ens_potext const v) {
  sim->ens.potext = v;
}

double ens_kindf(struct ens const* const ens) {
  return (double) ens->nmemb.dim / (2.0 * ens->epsilon);
}

double ens_kin_polybead_bw(struct ens const* const ens,
    size_t const ipoly, size_t const ibead) {
  double const m = ens->r[ipoly].mass / (2.0 * gsl_pow_2(ens->epsilon));

  if (ibead == 0) {
    size_t const jpoly = ens->r[ipoly].ifrom;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else {
      size_t const jbead = ens->nmemb.bead - 1;

      return m * ens_dist2(ens,
          &ens->r[ipoly].r[ibead], &ens->r[jpoly].r[jbead]);
    }
  } else {
    size_t const jbead = ibead - 1;

    return m * ens_dist2(ens,
        &ens->r[ipoly].r[ibead], &ens->r[ipoly].r[jbead]);
  }
}

double ens_kin_polybead_fw(struct ens const* const ens,
    size_t const ipoly, size_t const ibead) {
  double const m = ens->r[ipoly].mass / (2.0 * gsl_pow_2(ens->epsilon));

  if (ibead == ens->nmemb.bead - 1) {
    size_t const jpoly = ens->r[ipoly].ito;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else {
      size_t const jbead = 0;

      return m * ens_dist2(ens,
          &ens->r[ipoly].r[ibead], &ens->r[jpoly].r[jbead]);
    }
  } else {
    size_t const jbead = ibead + 1;

    return m * ens_dist2(ens,
        &ens->r[ipoly].r[ibead], &ens->r[ipoly].r[jbead]);
  }
}

double ens_kin_bead_fw(struct ens const* const ens, size_t const ibead) {
  double k = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    k += ens_kin_polybead_fw(ens, ipoly, ibead);

  return k;
}

double ens_kin_bead_bw(struct ens const* const ens, size_t const ibead) {
  double k = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    k += ens_kin_polybead_bw(ens, ipoly, ibead);

  return k;
}

double ens_kin_polybead(struct ens const* const ens,
    size_t const ipoly, size_t const ibead) {
  return ens_kin_polybead_bw(ens, ipoly, ibead) +
    ens_kin_polybead_fw(ens, ipoly, ibead);
}

double ens_kin_poly(struct ens const* const ens, size_t const ipoly) {
  double k = ens_kin_polybead_bw(ens, ipoly, 0);

  for (size_t ibead = 0; ibead < ens->nmemb.bead; ++ibead)
    k += ens_kin_polybead_fw(ens, ipoly, ibead);

  return k;
}

double ens_kin_total(struct ens const* const ens) {
  double k = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    k += ens_kin_poly(ens, ipoly);

  return k;
}

double ens_potint_bead(struct ens const* const ens, size_t const ibead) {
  double v = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < ens->nmemb.poly; ++jpoly)
      v += ens->potint(ens, &ens->r[ipoly].r[ibead], &ens->r[jpoly].r[ibead]);

  return v;
}

double ens_potext_polybead(struct ens const* const ens,
    size_t const ipoly, size_t const ibead) {
  return ens->potext(ens, &ens->r[ipoly].r[ibead]);
}

double ens_potext_bead(struct ens const* const ens, size_t const ibead) {
  double v = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    v += ens_potext_polybead(ens, ipoly, ibead);

  return v;
}

double ens_pot_bead(struct ens const* const ens, size_t const ibead) {
  return ens_potint_bead(ens, ibead) + ens_potext_bead(ens, ibead);
}

double ens_pot_total(struct ens const* const ens) {
  double v = 0.0;

  for (size_t ibead = 0; ibead < ens->nmemb.bead; ++ibead)
    v += ens_pot_bead(ens, ibead);

  return v;
}

double ens_est_pimc_td(struct ens const* const ens,
    __attribute__ ((__unused__)) void const* const p) {
  double const e = ens_kindf(ens);
  double const k = ens_kin_total(ens) / (double) ens->nmemb.bead;
  double const v = ens_pot_total(ens) / (double) ens->nmemb.bead;

  return e - k + v;
}

// TODO This is probably slightly wrong.
double ens_est_pigs_virial(struct ens const* const ens, void const* const p) {
  size_t const ibead = *(size_t const*) p;

  double const e = ens_kindf(ens);
  double k = 0.0;
  double const v = ens_pot_bead(ens, ibead);

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    k += ens_kin_polybead(ens, ipoly, ibead) / 2.0;

  return e - k + v;
}

// TODO This does not work with zero external potential.
double ens_est_pigs_mixed(struct ens const* const ens, void const* const p) {
  size_t const ibead = *(size_t const*) p;

  double v = 0.0;

  for (size_t ipoly = 0; ipoly < ens->nmemb.poly; ++ipoly)
    v += ens->potext(ens, &ens->r[ipoly].r[ibead]);

  return v;
}

void sim_weight_const(struct sim* const sim, void const* const p) {
  double const m = *(double const*) p;

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
    sim->ens.r[ipoly].mass = m;
}

void sim_perm_close(struct sim* const sim,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    sim->ens.r[ipoly].ifrom = ipoly;
    sim->ens.r[ipoly].ito = ipoly;
  }
}

void sim_perm_open(struct sim* const sim,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    sim->ens.r[ipoly].ifrom = SIZE_MAX;
    sim->ens.r[ipoly].ito = SIZE_MAX;
  }
}

void sim_perm_chain(struct sim* const sim,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    sim->ens.r[ipoly].ifrom = size_uwrap_dec(ipoly, sim->ens.nmemb.poly);
    sim->ens.r[ipoly].ito = size_uwrap_inc(ipoly, sim->ens.nmemb.poly);
  }
}

void sim_perm_random(struct sim* const sim,
    __attribute__ ((__unused__)) void const* const p) {
  size_t i[POLY_MAX];
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
    i[ipoly] = ipoly;

  gsl_ran_shuffle(sim->rng, i, sim->ens.nmemb.poly, sizeof *i);

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    sim->ens.r[i[ipoly]].ifrom = ipoly;
    sim->ens.r[ipoly].ito = i[ipoly];
  }
}

void sim_placer_point(struct sim* const sim, size_t const ipoly,
    struct bead const* const r, __attribute__ ((__unused__)) double const d,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      sim->ens.r[ipoly].r[ibead].r[idim] = r->r[idim];
}

void sim_placer_random(struct sim* const sim, size_t const ipoly,
    struct bead const* const r, double const d,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      sim->ens.r[ipoly].r[ibead].r[idim] = r->r[idim] +
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

      sim->ens.r[ipoly].r[ibead].r[idim] = r->r[idim] +
        d * cos(M_2PI * (k[idim] * t + phi)) / 2.0;
    }
  }
}

void sim_place_point(struct sim* const sim, sim_placer const f,
    void const* const p) {
  struct bead r;
  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
    r.r[idim] = 0.0;

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
    f(sim, ipoly, &r, sim->ens.length, p);
}

void sim_place_random(struct sim* const sim, sim_placer const f,
    void const* const p) {
  double const b = fp_rt((double) sim->ens.nmemb.poly,
      (double) sim->ens.nmemb.dim);

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    struct bead r;
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      r.r[idim] = ran_open(sim->rng, sim->ens.length / 2.0);

    f(sim, ipoly, &r, sim->ens.length / b, p);
  }
}

void sim_place_lattice(struct sim* const sim, sim_placer const f,
    void const* const p) {
  size_t const nlin = size_cirt(sim->ens.nmemb.poly, sim->ens.nmemb.dim);

  double const a = -sim->ens.length / 2.0;
  double const b = (double) nlin;

  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, nlin, i, sim->ens.nmemb.dim);

    struct bead r;
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      r.r[idim] = fp_lerp((double) i[idim], 0.0, b, -a, a);

    f(sim, ipoly, &r, sim->ens.length / b, p);
  }
}

void sim_move_null(__attribute__ ((__unused__)) struct sim* const sim) {}

void sim_move_accept_ss(struct sim* const sim) {
  ++sim->params.ssm.acc;
}

void sim_move_reject_ss(struct sim* const sim) {
  size_t const ipoly = sim->cache.ssm.ipoly;
  size_t const ibead = sim->cache.ssm.ibead;

  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
    sim->ens.r[ipoly].r[ibead].r[idim] = sim->cache.ssm.r.r[idim];

  ++sim->params.ssm.rej;
}

void sim_move_adjust_ss(struct sim* const sim) {
  sim->params.ssm.h = fp_clamp(sim->params.ssm.h *
      fp_cbalance((double) sim->params.ssm.acc, (double) sim->params.ssm.rej),
      DBL_EPSILON, sim->ens.length);
}

void sim_move_ss(struct sim* const sim,
    size_t const ipoly, size_t const ibead) {
  sim->accept = sim_move_accept_ss;
  sim->reject = sim_move_reject_ss;
  sim->adjust = sim_move_adjust_ss;

  sim->cache.ssm.ipoly = ipoly;
  sim->cache.ssm.ibead = ibead;

  sim_wrapper const f = sim->ens.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim) {
    sim->cache.ssm.r.r[idim] = sim->ens.r[ipoly].r[ibead].r[idim];

    sim->ens.r[ipoly].r[ibead].r[idim] =
      f(sim->ens.r[ipoly].r[ibead].r[idim] +
          sim->params.ssm.h * ran_open(sim->rng, 1.0), sim->ens.length);
  }
}

void sim_move_accept_cmd(struct sim* const sim) {
  ++sim->params.cmd.acc;
}

void sim_move_reject_cmd(struct sim* const sim) {
  size_t const ipoly = sim->cache.cmd.ipoly;

  for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      sim->ens.r[ipoly].r[ibead].r[idim] =
        sim->cache.cmd.r.r[ibead].r[idim];

  ++sim->params.cmd.rej;
}

void sim_move_adjust_cmd(struct sim* const sim) {
  sim->params.cmd.h = fp_clamp(sim->params.cmd.h *
      fp_cbalance((double) sim->params.cmd.acc, (double) sim->params.cmd.rej),
      DBL_EPSILON, sim->ens.length);
}

void sim_move_cmd(struct sim* const sim, size_t const ipoly) {
  sim->accept = sim_move_accept_cmd;
  sim->reject = sim_move_reject_cmd;
  sim->adjust = sim_move_adjust_cmd;

  sim->cache.cmd.ipoly = ipoly;

  sim_wrapper const f = sim->ens.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim) {
    double const x = sim->params.cmd.h *
      ran_open(sim->rng, 1.0);

    for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead) {
      sim->cache.cmd.r.r[ibead].r[idim] =
        sim->ens.r[ipoly].r[ibead].r[idim];

      sim->ens.r[ipoly].r[ibead].r[idim] =
        f(sim->ens.r[ipoly].r[ibead].r[idim] + x, sim->ens.length);
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

  double const V0 = ens_potint_bead(&sim->ens, ibead) +
    ens_potext_polybead(&sim->ens, ipoly, ibead);
  double const K0 = ens_kin_polybead(&sim->ens, ipoly, ibead);

  sim_move_ss(sim, ipoly, ibead);

  double const V1 = ens_potint_bead(&sim->ens, ibead) +
    ens_potext_polybead(&sim->ens, ipoly, ibead);
  double const K1 = ens_kin_polybead(&sim->ens, ipoly, ibead);

  return sim->ens.epsilon * (K1 + V1 - K0 - V0);
}

double sim_propose_cmd(struct sim* const sim) {
  size_t const ipoly = ran_index(sim->rng, sim->ens.nmemb.poly);

  double const V0 = ens_pot_total(&sim->ens);

  sim_move_cmd(sim, ipoly);

  double const V1 = ens_pot_total(&sim->ens);

  return sim->ens.epsilon * (V1 - V0);
}

void sim_decide_mq(struct sim* const sim, double const DeltaS) {
  if (DeltaS > 0.0 && exp(-DeltaS) < ran_uopen(sim->rng, 1.0))
    sim->reject(sim);
  else
    sim->accept(sim);

  sim->adjust(sim);
}

bool sim_res_close(__attribute__ ((__unused__)) struct sim const* const sim,
    FILE* const fp) {
  if (fclose(fp) == EOF)
    return false;

  return true;
}

FILE* sim_res_open(struct sim const* const sim, char const* const str) {
  // Sanity check.
  if (strchr(str, '.') != NULL || strchr(str, '/') != NULL)
    return NULL;

  char buf[BUFSIZ];
  int const k = snprintf(buf, sizeof buf, "%s/%s.data", sim->id, str);
  if (k < 0 || (size_t) k >= sizeof buf)
    return NULL;

  return fopen(buf, "w");
}

bool sim_res_print(struct sim const* const sim, char const* const str,
    sim_printer const f, void const* const p) {
  FILE* const fp = sim_res_open(sim, str);
  if (fp == NULL)
    return false;

  bool const q = f(sim, fp, p);

  if (!sim_res_close(sim, fp))
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
  if (fprintf(fp, "%g\n", sim->ens.length) < 0)
    return false;

  return true;
}

bool print_pots(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  size_t const npt = size_pow(sim->ens.nmemb.subdiv, sim->ens.nmemb.dim);

  struct bead r0;
  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
    r0.r[idim] = 0.0;

  for (size_t ipt = 0; ipt < npt; ++ipt) {
    size_t i[DIM_MAX];
    size_hc(ipt, sim->ens.nmemb.subdiv, i, sim->ens.nmemb.dim);

    struct bead r;
    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      r.r[idim] = fp_lerp((double) i[idim],
          0.0, (double) sim->ens.nmemb.subdiv, sim->ens.a, sim->ens.b);

    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      if (fprintf(fp, "%g ", r.r[idim]) < 0)
        return false;

    if (fprintf(fp, "%g %g\n",
          sim->ens.potext(&sim->ens, &r),
          sim->ens.potint(&sim->ens, &r0, &r)) < 0)
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

bool print_ndiv(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  if (fprintf(fp, "%zu\n", sim->ens.nmemb.div) < 0)
    return false;

  return true;
}

bool print_energy_bead(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
    if (fprintf(fp, "%g %g %g %g %g\n", (double) ibead * sim->ens.epsilon,
          stats_mean(sim->tdei[ibead]), stats_sem(sim->tdei[ibead]),
          stats_mean(sim->mixed[ibead]), stats_sem(sim->mixed[ibead])) < 0)
      return false;

  return true;
}

bool print_energy(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  size_t const ibead = sim->ens.nmemb.bead / 2;

  if (fprintf(fp, "%zu %g %g %g %g %g %g %g %g %g\n", sim->ens.istep.prod,
        ens_est_pimc_td(&sim->ens, NULL),
        stats_mean(sim->tde), stats_sem(sim->tde),
        ens_est_pigs_virial(&sim->ens, &ibead),
        stats_mean(sim->tdei[ibead]), stats_sem(sim->tdei[ibead]),
        ens_est_pigs_mixed(&sim->ens, &ibead),
        stats_mean(sim->mixed[ibead]), stats_sem(sim->mixed[ibead])) < 0)
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
  size_t const nbin = hist_nbin(sim->p);

  for (size_t ibin = 0; ibin < nbin; ++ibin) {
    struct bead r;
    hist_unbin(sim->p, r.r, ibin);

    for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
      if (fprintf(fp, "%g ", r.r[idim]) < 0)
        return false;

    if (fprintf(fp, "%g\n", hist_normhits(sim->p, ibin)) < 0)
      return false;
  }

  return true;
}

bool print_raddist(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  size_t const nbin = hist_nbin(sim->g);

  for (size_t ibin = 0; ibin < nbin; ++ibin) {
    double r0;
    hist_funbin(sim->g, &r0, ibin);

    double r1;
    hist_cunbin(sim->g, &r1, ibin);

    double const v = fp_ballvol(r1, sim->ens.nmemb.dim) -
      fp_ballvol(r0, sim->ens.nmemb.dim);

    double r;
    hist_unbin(sim->g, &r, ibin);

    if (fprintf(fp, "%g %g\n", r, hist_normhits(sim->g, ibin) * v) < 0)
      return false;
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool print_polys_bead(struct sim const* const sim, FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  if (fprintf(fp, "%g", (double) iindex * sim->ens.epsilon) < 0)
    return false;

  for (size_t idim = 0; idim < sim->ens.nmemb.dim; ++idim)
    if (fprintf(fp, " %g", sim->ens.r[ipoly].r[ibead].r[idim]) < 0)
      return false;

  if (fprintf(fp, "\n") < 0)
    return false;

  return true;
}

bool print_polys(struct sim const* const sim, FILE* const fp,
    __attribute__ ((__unused__)) void const* const p) {
  for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
    for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead)
      if (!print_polys_bead(sim, fp, ibead, ipoly, ibead))
        return false;

    if (sim->ens.r[ipoly].ito != SIZE_MAX)
      if (!print_polys_bead(sim, fp,
            sim->ens.nmemb.bead, sim->ens.r[ipoly].ito, 0))
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

  if (fprintf(fp,
        "(i_T + i_P) / (n_T + n_P) = (%zu + %zu) / (%zu + %zu) = %zu %%\n",
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
  size_t const ibead = sim->ens.nmemb.bead / 2;

  if (fprintf(fp, "E_T = %g +- %g (kappa = %g)\n"
        "E_L = %g +- %g (kappa = %g)\n"
        "E_M = %g +- %g (kappa = %g)\n",
        // ens_est_pimc_td
        stats_mean(sim->tde),
        100.0 * stats_sem(sim->tde), 100.0,
        // ens_est_pigs_virial
        stats_mean(sim->tdei[ibead]),
        100.0 * stats_sem(sim->tdei[ibead]), 100.0,
        // ens_est_pigs_mixed
        stats_mean(sim->mixed[ibead]),
        100.0 * stats_sem(sim->mixed[ibead]), 100.0) < 0)
    return false;

  return true;
}

bool sim_fini_fs(struct sim const* const sim) {
  // Sanity check.
  if (strncmp(sim->id, "run", 3) != 0)
    err_abort(strncmp);

  if (symlink(sim->id, "run-temp") == -1)
    err_abort(symlink);

  if (rename("run-temp", "run-latest") == -1)
    err_abort(rename);

  return true;
}

bool sim_init_fs(struct sim const* const sim) {
  // Sanity check.
  if (strncmp(sim->id, "run", 3) != 0)
    err_abort(strncmp);

  if (mkdir(sim->id, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    err_abort(mkdir);

  if (symlink(sim->id, "run-temp") == -1)
    err_abort(symlink);

  if (rename("run-temp", "run-current") == -1)
    err_abort(rename);

  return true;
}

bool sim_save_const(struct sim const* const sim) {
  return sim_res_print(sim, "periodic", print_periodic, NULL) &&
    sim_res_print(sim, "length", print_length, NULL) &&
    sim_res_print(sim, "pots", print_pots, NULL) &&
    sim_res_print(sim, "ndim", print_ndim, NULL) &&
    sim_res_print(sim, "npoly", print_npoly, NULL) &&
    sim_res_print(sim, "nbead", print_nbead, NULL) &&
    sim_res_print(sim, "nsubdiv", print_nsubdiv, NULL) &&
    sim_res_print(sim, "ndiv", print_ndiv, NULL);
}

bool sim_save_mut(struct sim const* const sim) {
  return sim_res_print(sim, "energy-bead", print_energy_bead, NULL) &&
    sim_res_print(sim, "energy-corrtime", print_energy_corrtime, NULL) &&
    sim_res_print(sim, "posdist", print_posdist, NULL) &&
    sim_res_print(sim, "raddist", print_raddist, NULL) &&
    sim_res_print(sim, "polys", print_polys, NULL);
}

void sim_free(struct sim* const sim) {
  if (sim != NULL) {
    hist_free(sim->g);
    hist_free(sim->p);

    for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead) {
      stats_free(sim->mixed[ibead]);
      stats_free(sim->tdei[ibead]);
    }
    stats_free(sim->tde);

    gsl_rng_free(sim->rng);
  }

  free(sim);
}

struct sim* sim_alloc(size_t const ndim, size_t const npoly,
    size_t const nbead, size_t const nsubdiv, size_t const ndiv,
    size_t const nthrm, size_t const nprod,
    size_t const nthrmrec, size_t const nprodrec,
    bool const periodic,
    double const length, double const mass, double const betau) {
  dynamic_assert(ndim <= DIM_MAX, "too many dimensions");
  dynamic_assert(npoly <= POLY_MAX, "too many polymers");
  dynamic_assert(nbead <= BEAD_MAX, "too many beads");
  dynamic_assert(nsubdiv <= SUBDIV_MAX, "too many subdivisions");
  dynamic_assert(ndiv <= DIV_MAX, "too many divisions");

  // These assertions guarantee that adding or multiplying two steps or indices
  // never wraps around, which makes their manipulation easier.
  dynamic_assert(nthrm <= SQRT_SIZE_MAX, "too many thermalization steps");
  dynamic_assert(nprod <= SQRT_SIZE_MAX, "too many production steps");
  dynamic_assert(nthrmrec <= nthrm, "too many thermalization recording steps");
  dynamic_assert(nprodrec <= nprod, "too many production recording steps");

  dynamic_assert(length > 0.0, "empty simulation box");
  dynamic_assert(mass > 0.0, "negative default mass");
  dynamic_assert(betau > 0.0, "weird things going on");

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

    sim->ens.periodic = periodic;
    sim->ens.length = length;

    if (periodic) {
      sim->ens.a = 0.0;
      sim->ens.b = length;
    } else {
      sim->ens.a = -length / 2.0;
      sim->ens.b = length / 2.0;
    }

    sim->ens.epsilon = betau / (double) nbead;

    sim->ens.potint = ens_pot_zero;
    sim->ens.potext = ens_potext_zero;

    sim->ens.nmemb.dim = ndim;
    sim->ens.nmemb.poly = npoly;
    sim->ens.nmemb.bead = nbead;
    sim->ens.nmemb.subdiv = nsubdiv;
    sim->ens.nmemb.div = ndiv;

    sim->ens.nstep.thrm = nthrm;
    sim->ens.nstep.prod = nprod;
    sim->ens.nstep.thrmrec = nthrmrec;
    sim->ens.nstep.prodrec = nprodrec;

    sim->ens.istep.thrm = 0;
    sim->ens.istep.prod = 0;
    sim->ens.istep.thrmrec = 0;
    sim->ens.istep.prodrec = 0;

    sim->accept = sim_move_null;
    sim->reject = sim_move_null;
    sim->adjust = sim_move_null;

    sim->tde = stats_alloc(nprod, true);
    if (sim->tde == NULL)
      p = false;

    for (size_t ibead = 0; ibead < nbead; ++ibead) {
      sim->tdei[ibead] = stats_alloc(nprod, false);
      if (sim->tdei[ibead] == NULL)
        p = false;

      sim->mixed[ibead] = stats_alloc(nprod, false);
      if (sim->mixed[ibead] == NULL)
        p = false;
    }

    sim->p = hist_alloc(ndim, nsubdiv, sim->ens.a, sim->ens.b);
    if (sim->p == NULL)
      p = false;

    sim->g = hist_alloc(1, ndiv, 0.0, length / 2.0);
    if (sim->g == NULL)
      p = false;

    sim->params.ssm.acc = 0;
    sim->params.ssm.rej = 0;
    sim->params.ssm.h = length / 2.0;

    sim->params.cmd.acc = 0;
    sim->params.cmd.rej = 0;
    sim->params.cmd.h = length / 2.0;

    sim_weight_const(sim, &mass);
    sim_perm_close(sim, NULL);
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

  if (!sim_init_fs(sim))
    err_abort(sim_init_fs);

  if (!sim_save_const(sim))
    err_abort(sim_save_const);

  printf("Starting simulation '%s'...\n", sim->id);

  FILE* const paramsfp = sim_res_open(sim, "params");
  if (paramsfp == NULL)
    err_abort(sim_res_open);

  FILE* const energyfp = sim_res_open(sim, "energy");
  if (energyfp == NULL)
    err_abort(sim_res_open);

  sim_proposer const proposers[] = {sim_propose_ss, sim_propose_cmd};
  size_t const nproposer = sizeof proposers / sizeof *proposers;

  for (size_t istep = 0;
      istep < sim->ens.nstep.thrm + sim->ens.nstep.prod;
      ++istep) {
    sim_decide_mq(sim, proposers[ran_index(sim->rng, nproposer)](sim));

    int signum;
    if (sigs_use(&signum))
      switch (signum) {
        case SIGUSR1:
          (void) print_progress(sim, stdout, NULL);
          break;
        case SIGUSR2:
          (void) sim_save_mut(sim);
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
      (void) stats_accum(sim->tde, ens_est_pimc_td(&sim->ens, NULL));

      for (size_t ibead = 0; ibead < sim->ens.nmemb.bead; ++ibead) {
        (void) stats_accum(sim->tdei[ibead],
            ens_est_pigs_virial(&sim->ens, &ibead) / sim->ens.nmemb.poly);
        (void) stats_accum(sim->mixed[ibead],
            ens_est_pigs_mixed(&sim->ens, &ibead) / sim->ens.nmemb.poly);
      }

      for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly) {
        size_t const ibead = sim->ens.nmemb.bead / 2;

        (void) hist_accum(sim->p, sim->ens.r[ipoly].r[ibead].r);
      }

      for (size_t ipoly = 0; ipoly < sim->ens.nmemb.poly; ++ipoly)
        for (size_t jpoly = ipoly + 1; jpoly < sim->ens.nmemb.poly; ++jpoly) {
          size_t const ibead = sim->ens.nmemb.bead / 2;

          double const d = ens_dist(&sim->ens,
              &sim->ens.r[ipoly].r[ibead], &sim->ens.r[jpoly].r[ibead]);

          (void) hist_accum(sim->g, &d);
        }

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

  if (!sim_res_close(sim, energyfp))
    err_abort(sim_res_close);

  if (!sim_res_close(sim, paramsfp))
    err_abort(sim_res_close);

  printf("Saving results of simulation '%s'...\n", sim->id);

  if (!sim_save_mut(sim))
    err_abort(sim_save_mut);

  if (!print_progress(sim, stdout, NULL))
    err_abort(print_progress);

  if (!print_wrong_results_fast(sim, stdout, NULL))
    err_abort(print_wrong_results_fast);

  if (!sim_fini_fs(sim))
    err_abort(sim_fini_fs);

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    err_abort(feenableexcept);
#endif
#endif

  return true;
}
