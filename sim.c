#include "err.h"
#include "exts.h"
#include "fp.h"
#include "lims.h"
#include "ran.h"
#include "sigs.h"
#include "size.h"
#include "stats.h"
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

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

// This structure represents a bead and is a part of `poly`.
struct bead {
  double d[DIM_MAX];
};

// This structure represents a polymer and is a part of `ensem`.
struct poly {
  double m;
  size_t ifrom;
  size_t ito;
  struct bead r[BEAD_MAX];
};

// This structure represents an ensemble and is a part of `napkin`.
struct ensem {
  bool periodic;
  double L;
  double beta;
  double tau;
  double (* Vint)(struct ensem const*, struct bead const*, struct bead const*);
  double (* Vend)(struct ensem const*, struct bead const*, struct bead const*);
  double (* Vext)(struct ensem const*, struct bead const*);
  struct memb Nmemb;
  struct step Nstep;
  struct step istep;
  struct poly R[POLY_MAX];
};

// This structure contains all the necessary bookkeeping.
struct napkin {
  gsl_rng* rng;
  struct stats* tde;
  size_t* p;
  size_t* g;
  void (* accept)(struct napkin*);
  void (* reject)(struct napkin*);
  void (* adjust)(struct napkin*);
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
      size_t accepted;
      size_t rejected;
      double dx;
    } ssm;
    struct {
      size_t accepted;
      size_t rejected;
      double dx;
    } cmd;
  } params;
  struct ensem ensem;
  char id[BUFSIZ];
};

double bead_norm2(struct ensem const* const ensem,
    struct bead const* const r) {
  double s = 0.0;

  double (* const f)(double, double) =
    ensem->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ensem->Nmemb.dim; ++idim)
    s += gsl_pow_2(f(r->d[idim], ensem->L));

  return s;
}

double bead_norm(struct ensem const* const ensem,
    struct bead const* const r) {
  return sqrt(bead_norm2(ensem, r));
}

double bead_dist2(struct ensem const* const ensem,
    struct bead const* const r0, struct bead const* const r1) {
  double s = 0.0;

  double (* const f)(double, double) =
    ensem->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ensem->Nmemb.dim; ++idim)
    s += gsl_pow_2(f(r1->d[idim] - r0->d[idim], ensem->L));

  return s;
}

double bead_dist(struct ensem const* const ensem,
    struct bead const* const r0, struct bead const* const r1) {
  return sqrt(bead_dist2(ensem, r0, r1));
}

// Kinetic and Potential Energies (K, V) --------------------------------------

// Kinetic energy within a polymer from one bead to the previous one.
__attribute__ ((__nonnull__, __pure__))
static double K_polybead_bw(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  double const M = ensem->R[ipoly].m / (2.0 * gsl_pow_2(ensem->tau));

  if (ibead == 0) {
    size_t const jpoly = ensem->R[ipoly].ifrom;
    size_t const jbead = ensem->Nmemb.bead - 1;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return M * bead_dist2(ensem,
          &ensem->R[ipoly].r[ibead], &ensem->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead - 1;

    return M * bead_dist2(ensem,
        &ensem->R[ipoly].r[ibead], &ensem->R[ipoly].r[jbead]);
  }
}

// Kinetic energy within a polymer from one bead to the next one.
// \lambda = \hbar / 2 m, K = x^2 / 4 \lambda \tau^2 = M x^2
// M = m / 2 \hbar \tau^2
__attribute__ ((__nonnull__, __pure__))
static double K_polybead_fw(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  double const M = ensem->R[ipoly].m / (2.0 * gsl_pow_2(ensem->tau));

  if (ibead == ensem->Nmemb.bead - 1) {
    size_t const jpoly = ensem->R[ipoly].ito;
    size_t const jbead = 0;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return M * bead_dist2(ensem,
          &ensem->R[ipoly].r[ibead], &ensem->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead + 1;

    return M * bead_dist2(ensem,
        &ensem->R[ipoly].r[ibead], &ensem->R[ipoly].r[jbead]);
  }
}

// Kinetic energy within a polymer from one bead to its two neighbors.
__attribute__ ((__nonnull__, __pure__))
static double K_polybead(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  return K_polybead_bw(ensem, ipoly, ibead) +
    K_polybead_fw(ensem, ipoly, ibead);
}

// Kinetic energy within a polymer.
__attribute__ ((__nonnull__, __pure__))
static double K_poly(struct ensem const* const ensem,
    size_t const ipoly) {
  double K = 0.0;

  for (size_t ibead = 0; ibead < ensem->Nmemb.bead; ++ibead)
    K += K_polybead_fw(ensem, ipoly, ibead);

  return K;
}

// Kinetic energy.
__attribute__ ((__nonnull__, __pure__))
static double K_total(struct ensem const* const ensem) {
  double K = 0.0;

  for (size_t ipoly = 0; ipoly < ensem->Nmemb.poly; ++ipoly)
    K += K_poly(ensem, ipoly);

  return K;
}

// Potential energy between certain beads of all polymers.
__attribute__ ((__nonnull__, __pure__))
static double Vint_bead(struct ensem const* const ensem,
    size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < ensem->Nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < ensem->Nmemb.poly; ++jpoly)
      V += ensem->Vint(ensem,
          &ensem->R[ipoly].r[ibead], &ensem->R[jpoly].r[ibead]);

  return V;
}

// Additional potential energy between either of the ends of all polymers.
__attribute__ ((__nonnull__, __pure__))
static double Vend_bead(struct ensem const* const ensem,
    size_t const ibead) {
  double V = 0.0;

  if (ibead == 0) {
    for (size_t ipoly = 0; ipoly < ensem->Nmemb.poly; ++ipoly)
      if (ensem->R[ipoly].ifrom == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < ensem->Nmemb.poly; ++jpoly)
          if (ensem->R[jpoly].ifrom == SIZE_MAX)
            V += ensem->Vend(ensem,
                &ensem->R[ipoly].r[ibead], &ensem->R[jpoly].r[ibead]);
  } else if (ibead == ensem->Nmemb.bead - 1)
    for (size_t ipoly = 0; ipoly < ensem->Nmemb.poly; ++ipoly)
      if (ensem->R[ipoly].ito == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < ensem->Nmemb.poly; ++jpoly)
          if (ensem->R[jpoly].ito == SIZE_MAX)
            V += ensem->Vend(ensem,
                &ensem->R[ipoly].r[ibead], &ensem->R[jpoly].r[ibead]);

  return V;
}

// External potential energy for one bead in a polymer.
__attribute__ ((__nonnull__, __pure__))
static double Vext_polybead(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  return ensem->Vext(ensem, &ensem->R[ipoly].r[ibead]);
}

// External potential energy for every bead in a polymer.
__attribute__ ((__nonnull__, __pure__))
static double Vext_bead(struct ensem const* const ensem, size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < ensem->Nmemb.poly; ++ipoly)
    V += Vext_polybead(ensem, ipoly, ibead);

  return V;
}

// Potential energy for one bead over all polymers.
__attribute__ ((__nonnull__, __pure__))
static double V_bead(struct ensem const* const ensem,
    size_t const ibead) {
  return Vint_bead(ensem, ibead) + Vend_bead(ensem, ibead) +
    Vext_bead(ensem, ibead);
}

// Potential energy.
__attribute__ ((__nonnull__, __pure__))
static double V_total(struct ensem const* const ensem) {
  double V = 0.0;

  for (size_t ibead = 0; ibead < ensem->Nmemb.bead; ++ibead)
    V += V_bead(ensem, ibead);

  return V;
}

// Estimators (est) -----------------------------------------------------------

// Thermodynamic energy estimator for closed polymers.
// \langle E_T\rangle & = \frac 1{M \tau} \sum_{k = 1}^M
// \Bigl\langle\frac{d N} 2 - \frac{|R_{k + 1 \bmod M} - R_k|^2}{4 \lambda \tau} + \tau V(R_k)\Bigr\rangle
__attribute__ ((__nonnull__, __pure__))
static double est_pimc_tde(struct ensem const* const ensem) {
  double const E = (double) (ensem->Nmemb.dim * ensem->Nmemb.poly) /
    (2.0 * ensem->tau);
  double const K = K_total(ensem) /
    (double) (ensem->Nmemb.poly * ensem->Nmemb.bead);
  double const V = V_total(ensem) /
    (double) (ensem->Nmemb.poly * ensem->Nmemb.bead);

  return E - K + V;
}

// TODO The crap thermodynamic energy estimator for open polymers.
__attribute__ ((__nonnull__, __pure__))
static double est_pigs_crap(struct ensem const* const ensem) {
  double const V = Vext_bead(ensem, ensem->Nmemb.bead / 2) /
    (double) ensem->Nmemb.poly;

  return NAN;
}

// TODO Thermodynamic energy estimator for open polymers.
__attribute__ ((__nonnull__, __pure__))
static double est_pigs_tde(struct ensem const* const ensem) {
  return NAN;
}

// Mass Configurations (Rm) ---------------------------------------------------

// Equal mass for every particle.
__attribute__ ((__nonnull__))
static void Rm_const(struct napkin* const napkin, double const m) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly)
    napkin->ensem.R[ipoly].m = m;
}

// End Configurations (Rend) --------------------------------------------------

// Close the current configuration.
__attribute__ ((__nonnull__))
static void Rend_close(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].ifrom = ipoly;
    napkin->ensem.R[ipoly].ito = ipoly;
  }
}

// Open the current configuration.
__attribute__ ((__nonnull__))
static void Rend_open(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].ifrom = SIZE_MAX;
    napkin->ensem.R[ipoly].ito = SIZE_MAX;
  }
}

// Cycle the current configuration.
__attribute__ ((__nonnull__))
static void Rend_cycle(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].ifrom = size_uwrap_dec(ipoly, napkin->ensem.Nmemb.poly);
    napkin->ensem.R[ipoly].ito = size_uwrap_inc(ipoly, napkin->ensem.Nmemb.poly);
  }
}

// Bead Configurations (Rr) ---------------------------------------------------

// Random initial configuration.
__attribute__ ((__nonnull__))
static void Rr_rand(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = ran_uopen(napkin->rng,
            napkin->ensem.L);
}

// Random point initial configuration.
__attribute__ ((__nonnull__))
static void Rr_randpt(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly)
    for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim) {
      double const x = ran_uopen(napkin->rng, napkin->ensem.L);

      for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = x;
    }
}

// Random lattice initial configuration.
__attribute__ ((__nonnull__))
static void Rr_randlatt(struct napkin* const napkin) {
  size_t const Nlin = size_cirt(napkin->ensem.Nmemb.poly,
      napkin->ensem.Nmemb.dim);

  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, Nlin, i, napkin->ensem.Nmemb.dim);

    for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = fp_lerp((double) i[idim] +
            ran_uopen(napkin->rng, 1.0),
            0.0, (double) Nlin,
            0.0, napkin->ensem.L);
  }
}

// Point lattice initial configuration.
__attribute__ ((__nonnull__))
static void Rr_ptlatt(struct napkin* const napkin) {
  size_t const Nlin = size_cirt(napkin->ensem.Nmemb.poly,
      napkin->ensem.Nmemb.dim);

  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, Nlin, i, napkin->ensem.Nmemb.dim);

    for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = fp_lerp((double) i[idim],
            0.0, (double) Nlin,
            0.0, napkin->ensem.L);
  }
}

// Circle lattice initial configuration.
// The equations are based on the theory of Lissajous knots.
__attribute__ ((__nonnull__))
static void Rr_circlatt(struct napkin* const napkin) {
  size_t const Nlin = size_cirt(napkin->ensem.Nmemb.poly,
      napkin->ensem.Nmemb.dim);

  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, Nlin, i, napkin->ensem.Nmemb.dim);

    for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead) {
      double const t = (double) ibead / (double) napkin->ensem.Nmemb.bead;

      for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim) {
        double const phi = (double) idim / (double) napkin->ensem.Nmemb.dim;
        double const x = sin(M_2PI * (t + phi / 2.0));

        napkin->ensem.R[ipoly].r[ibead].d[idim] = fp_lerp((double) i[idim] +
            0.5 + x / 4.0,
            0.0, (double) Nlin,
            0.0, napkin->ensem.L);
      }
    }
  }
}

// Physical Moves (move) ------------------------------------------------------

__attribute__ ((__nonnull__))
static void move_accept_ss(struct napkin* const napkin) {
  ++napkin->params.ssm.accepted;
}

__attribute__ ((__nonnull__))
static void move_reject_ss(struct napkin* const napkin) {
  size_t const ipoly = napkin->hist.ssm.ipoly;
  size_t const ibead = napkin->hist.ssm.ibead;

  for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
    napkin->ensem.R[ipoly].r[ibead].d[idim] = napkin->hist.ssm.r.d[idim];

  ++napkin->params.ssm.rejected;
}

__attribute__ ((__nonnull__))
static void move_adjust_ss(struct napkin* const napkin) {
  napkin->params.ssm.dx = napkin->params.ssm.dx *
    fp_dbalance((double) napkin->params.ssm.accepted,
        (double) napkin->params.ssm.rejected);
}

__attribute__ ((__nonnull__))
static void move_ss(struct napkin* const napkin,
    size_t const ipoly, size_t const ibead) {
  napkin->accept = move_accept_ss;
  napkin->reject = move_reject_ss;
  napkin->adjust = move_adjust_ss;

  napkin->hist.ssm.ipoly = ipoly;
  napkin->hist.ssm.ibead = ibead;

  double (* const f)(double, double) =
    napkin->ensem.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim) {
    napkin->hist.ssm.r.d[idim] = napkin->ensem.R[ipoly].r[ibead].d[idim];

    napkin->ensem.R[ipoly].r[ibead].d[idim] =
      f(napkin->ensem.R[ipoly].r[ibead].d[idim] +
        napkin->params.ssm.dx * ran_open(napkin->rng, 1.0), napkin->ensem.L);
  }
}

__attribute__ ((__nonnull__))
static void move_accept_cmd(struct napkin* const napkin) {
  ++napkin->params.cmd.accepted;
}

__attribute__ ((__nonnull__))
static void move_reject_cmd(struct napkin* const napkin) {
  size_t const ipoly = napkin->hist.cmd.ipoly;

  for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
      napkin->ensem.R[ipoly].r[ibead].d[idim] =
        napkin->hist.cmd.R.r[ibead].d[idim];

  ++napkin->params.cmd.rejected;
}

__attribute__ ((__nonnull__))
static void move_adjust_cmd(struct napkin* const napkin) {
  napkin->params.cmd.dx = napkin->params.cmd.dx *
    fp_cbalance((double) napkin->params.cmd.accepted,
        (double) napkin->params.cmd.rejected);
}

__attribute__ ((__nonnull__))
static void move_cmd(struct napkin* const napkin,
    size_t const ipoly) {
  napkin->accept = move_accept_cmd;
  napkin->reject = move_reject_cmd;
  napkin->adjust = move_adjust_cmd;

  napkin->hist.cmd.ipoly = ipoly;

  double (* const f)(double, double) =
    napkin->ensem.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim) {
    double const x = napkin->params.cmd.dx *
      ran_open(napkin->rng, 1.0);

    for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead) {
      napkin->hist.cmd.R.r[ibead].d[idim] =
        napkin->ensem.R[ipoly].r[ibead].d[idim];

      napkin->ensem.R[ipoly].r[ibead].d[idim] =
        f(napkin->ensem.R[ipoly].r[ibead].d[idim] + x, napkin->ensem.L);
    }
  }
}

__attribute__ ((__noreturn__))
static void move_accept_bisect(
    __attribute__ ((__unused__)) struct napkin* const napkin) {
  err_abort(NULL);
}

__attribute__ ((__noreturn__))
static void move_reject_bisect(
    __attribute__ ((__unused__)) struct napkin* const napkin) {
  err_abort(NULL);
}

__attribute__ ((__noreturn__))
static void move_bisect(
    __attribute__ ((__unused__)) struct napkin* const napkin,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

__attribute__ ((__noreturn__))
static void move_accept_swap(
    __attribute__ ((__unused__)) struct napkin* const napkin) {
  err_abort(NULL);
}

__attribute__ ((__noreturn__))
static void move_reject_swap(
    __attribute__ ((__unused__)) struct napkin* const napkin) {
  err_abort(NULL);
}

__attribute__ ((__noreturn__))
static void move_swap(
    __attribute__ ((__unused__)) struct napkin* const napkin,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

// Workers (work) -------------------------------------------------------------

__attribute__ ((__nonnull__))
static double work_ss(struct napkin* const napkin) {
  size_t const ipoly = ran_index(napkin->rng, napkin->ensem.Nmemb.poly);
  size_t const ibead = ran_index(napkin->rng, napkin->ensem.Nmemb.bead);

  double const V0 =
    Vint_bead(&napkin->ensem, ibead) + Vend_bead(&napkin->ensem, ibead) +
    Vext_polybead(&napkin->ensem, ipoly, ibead);
  double const K0 = K_polybead(&napkin->ensem, ipoly, ibead);

  move_ss(napkin, ipoly, ibead);

  double const V1 =
    Vint_bead(&napkin->ensem, ibead) + Vend_bead(&napkin->ensem, ibead) +
    Vext_polybead(&napkin->ensem, ipoly, ibead);
  double const K1 = K_polybead(&napkin->ensem, ipoly, ibead);

  return napkin->ensem.tau * (K1 + V1 - K0 - V0);
}

__attribute__ ((__nonnull__))
static double work_cmd(struct napkin* const napkin) {
  size_t const ipoly = ran_index(napkin->rng, napkin->ensem.Nmemb.poly);

  double const V0 = V_total(&napkin->ensem);

  move_cmd(napkin, ipoly);

  double const V1 = V_total(&napkin->ensem);

  return napkin->ensem.tau * (V1 - V0);
}

// Work Horses () -------------------------------------------------------------

__attribute__ ((__nonnull__))
static void choose(struct napkin* const napkin,
    double const DeltaS) {
  if (DeltaS <= 0.0 || ran_uopen(napkin->rng, 1.0) <= exp(-DeltaS))
    napkin->accept(napkin);
  else
    napkin->reject(napkin);

  napkin->adjust(napkin);
}

__attribute__ ((__nonnull__))
static void prob_accum(struct napkin const* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead) {
      size_t i[DIM_MAX];

      for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
        i[idim] = size_uclamp(
            (size_t) (floor(fp_lerp(napkin->ensem.R[ipoly].r[ibead].d[idim],
                  0.0, (double) napkin->ensem.L,
                  0.0, (double) napkin->ensem.Nmemb.subdiv))),
            napkin->ensem.Nmemb.subdiv);

      ++napkin->p[size_unhc(napkin->ensem.Nmemb.subdiv, i,
          napkin->ensem.Nmemb.dim)];
    }
}

__attribute__ ((__nonnull__))
static void paircorr_accum(struct napkin const* const napkin) {
  size_t const Nbin = size_pow(napkin->ensem.Nmemb.subdiv,
      napkin->ensem.Nmemb.dim);

  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < napkin->ensem.Nmemb.poly; ++jpoly)
      for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead) {
        size_t const i = size_uclamp(
            (size_t) (floor(fp_lerp(bead_dist(&napkin->ensem,
                    &napkin->ensem.R[ipoly].r[ibead],
                    &napkin->ensem.R[jpoly].r[ibead]),
                  0.0, (double) napkin->ensem.L,
                  0.0, (double) Nbin))), Nbin);

        ++napkin->g[i];
      }
}

// Good New Printing Stuff (disp) ---------------------------------------------

__attribute__ ((__nonnull__))
static bool res_close(
    __attribute__ ((__unused__)) struct napkin const* const napkin,
    FILE* const fp) {
  if (fclose(fp) == EOF)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static FILE* res_open(struct napkin const* const napkin,
    char const* const str) {
  // Sanity check.
  if (strchr(str, '.') != NULL || strchr(str, '/') != NULL)
    return NULL;

  char buf[BUFSIZ];
  int const k = snprintf(buf, sizeof buf,
      "%s/%s.data", napkin->id, str);
  if (k < 0 || (size_t) k >= sizeof buf)
    return NULL;

  return fopen(buf, "w");
}

__attribute__ ((__nonnull__))
static bool res_use(struct napkin const* const napkin, char const* const str,
    bool (* const f)(struct napkin const*, FILE*)) {
  FILE* const fp = res_open(napkin, str);
  if (fp == NULL)
    return false;

  bool const p = f(napkin, fp);

  if (!res_close(napkin, fp))
    return false;

  return p;
}

__attribute__ ((__nonnull__))
static bool disp_ndim(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.Nmemb.dim) < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_npoly(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.Nmemb.poly) < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_nbead(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.Nmemb.bead) < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_nsubdiv(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.Nmemb.subdiv) < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_L(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%f\n", napkin->ensem.L) < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_R_r(struct napkin const* const napkin, FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  if (fprintf(fp, "%zu", iindex) < 0)
    return false;

  for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
    if (fprintf(fp, " %f", napkin->ensem.R[ipoly].r[ibead].d[idim]) < 0)
      return false;

  if (fprintf(fp, "\n") < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_R(struct napkin const* const napkin, FILE* const fp) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.Nmemb.poly; ++ipoly) {
    for (size_t ibead = 0; ibead < napkin->ensem.Nmemb.bead; ++ibead)
      if (!disp_R_r(napkin, fp, ibead, ipoly, ibead))
        return false;

    if (napkin->ensem.R[ipoly].ito != SIZE_MAX)
      if (!disp_R_r(napkin, fp,
            napkin->ensem.Nmemb.bead, napkin->ensem.R[ipoly].ito, 0))
        return false;

    if (fprintf(fp, "\n") < 0)
      return false;
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_params(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu %zu %zu %zu %f %zu %zu %f\n",
        napkin->ensem.istep.thrm, napkin->ensem.istep.prod,
        napkin->params.ssm.accepted, napkin->params.ssm.rejected,
        napkin->params.ssm.dx,
        napkin->params.cmd.accepted, napkin->params.cmd.rejected,
        napkin->params.cmd.dx) < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_E(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu %f %f %f\n",
        napkin->ensem.istep.prod, est_pimc_tde(&napkin->ensem),
        stats_mean(napkin->tde), stats_sem(napkin->tde)) < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_V(struct napkin const* const napkin, FILE* const fp) {
  size_t const Nbin = size_pow(napkin->ensem.Nmemb.subdiv,
      napkin->ensem.Nmemb.dim);

  size_t const Npt = size_pow(napkin->ensem.Nmemb.subdiv + 1,
      napkin->ensem.Nmemb.dim);

  struct bead r0;
  for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
    r0.d[idim] = 0;

  for (size_t ipt = 0; ipt < Npt; ++ipt) {
    size_t i[DIM_MAX];
    size_hc(ipt, napkin->ensem.Nmemb.subdiv + 1, i, napkin->ensem.Nmemb.dim);

    struct bead r;

    for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim],
          0.0, (double) napkin->ensem.Nmemb.subdiv,
          0.0, napkin->ensem.L);

    for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
      if (fprintf(fp, "%f ", r.d[idim]) < 0)
        return false;

    // TODO Unclamp.
    if (fprintf(fp, "%f %f\n",
          napkin->ensem.Vext(&napkin->ensem, &r),
          fp_clamp(napkin->ensem.Vint(&napkin->ensem, &r0, &r),
            -8.0, 8.0)) < 0)
      return false;
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_p(struct napkin const* const napkin, FILE* const fp) {
  size_t const Nbin = size_pow(napkin->ensem.Nmemb.subdiv,
      napkin->ensem.Nmemb.dim);

  size_t N = 0;
  for (size_t ibin = 0; ibin < Nbin; ++ibin)
    N += napkin->p[ibin];
  double const S = (double) N * napkin->ensem.L;

  for (size_t ibin = 0; ibin < Nbin; ++ibin) {
    size_t i[DIM_MAX];
    size_hc(ibin, napkin->ensem.Nmemb.subdiv, i, napkin->ensem.Nmemb.dim);

    struct bead r;

    for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim] + 0.5,
          0.0, (double) napkin->ensem.Nmemb.subdiv,
          0.0, napkin->ensem.L);

    for (size_t idim = 0; idim < napkin->ensem.Nmemb.dim; ++idim)
      if (fprintf(fp, "%f ", r.d[idim]) < 0)
        return false;

    if (fprintf(fp, "%f\n", (double) napkin->p[ibin] / S) < 0)
      return false;
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_g(struct napkin const* const napkin, FILE* const fp) {
  size_t const Nbin = size_pow(napkin->ensem.Nmemb.subdiv,
      napkin->ensem.Nmemb.dim);

  size_t N = 0.0;
  for (size_t ibin = 0; ibin < Nbin; ++ibin)
    N += napkin->g[ibin];
  double const S = N == 0 ? 1.0 : (double) N;

  for (size_t ibin = 0; ibin < Nbin; ++ibin) {
    double const r = fp_lerp((double) ibin + 0.5,
        0.0, (double) napkin->ensem.Nmemb.subdiv,
        0.0, napkin->ensem.L);

    if (fprintf(fp, "%f %f\n", r, (double) napkin->g[ibin] / S) < 0)
      return false;
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool disp_i(struct napkin const* const napkin, FILE* const fp) {
  if (printf("i / N = (%zu + %zu) / (%zu + %zu) = %zu %%\n",
      napkin->ensem.istep.thrm, napkin->ensem.istep.prod, napkin->ensem.Nstep.thrm, napkin->ensem.Nstep.prod,
      100 * (napkin->ensem.istep.thrm + napkin->ensem.istep.prod) /
      (napkin->ensem.Nstep.thrm + napkin->ensem.Nstep.prod)) < 0)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool gen_id(gsl_rng* const rng, char* const buf, size_t const n) {
  time_t t;
  if (time(&t) == (time_t) -1)
    return false;

  struct tm tm;
  if (localtime_r(&t, &tm) == NULL)
    return false;

  int const k = snprintf(buf, n,
      "run-%04d-%02d-%02d-%02d-%02d-%02d-%03zx",
      tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
      tm.tm_hour, tm.tm_min, tm.tm_sec,
      ran_index(rng, 0x1000));
  if (k < 0 || (size_t) k >= n)
    return false;

  return true;
}

__attribute__ ((__nonnull__))
static bool prepare_stuff(struct napkin const* const napkin) {
  // Sanity check.
  if (strncmp(napkin->id, "run", 3) != 0)
    err_abort(strncmp);

  if (symlink(napkin->id, napkin->id) == -1)
    err_abort(symlink);

  if (rename(napkin->id, "run-latest") == -1)
    err_abort(rename);

  if (mkdir(napkin->id, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    err_abort(mkdir);

  return true;
}

// Other Crap () --------------------------------------------------------------

__attribute__ ((__nonnull__))
static void show_progress(struct napkin const* const napkin) {
  (void) printf("E = %f +- %f (kappa = 0)\n",
      napkin->ensem.istep.thrm, napkin->ensem.istep.prod, napkin->ensem.Nstep.thrm, napkin->ensem.Nstep.prod,
      100 * (napkin->ensem.istep.thrm + napkin->ensem.istep.prod) /
      (napkin->ensem.Nstep.thrm + napkin->ensem.Nstep.prod),
      stats_mean(napkin->tde), stats_sem(napkin->tde));
}

__attribute__ ((__nonnull__))
static void show_results(struct napkin const* const napkin) {
  (void) printf("E = %f +- %f (kappa = %f)\n",
      stats_mean(napkin->tde), stats_corrsem(napkin->tde),
      stats_corrtime(napkin->tde));
}

__attribute__ ((__nonnull__))
static void save_const(struct napkin const* const napkin) {
  (void) res_use(napkin, "length", disp_L);
  (void) res_use(napkin, "ndim", disp_ndim);
  (void) res_use(napkin, "npoly", disp_npoly);
  (void) res_use(napkin, "nbead", disp_nbead);
  (void) res_use(napkin, "nsubdiv", disp_nsubdiv);
  (void) res_use(napkin, "pots", disp_V);
}

__attribute__ ((__nonnull__))
static void save_mut(struct napkin const* const napkin) {
  (void) res_use(napkin, "polys", disp_R);
  (void) res_use(napkin, "paircorr", disp_g);
  (void) res_use(napkin, "posdist", disp_p);
}

// Napkin Management (napkin) -------------------------------------------------

static void napkin_free(struct napkin* const napkin) {
  if (napkin != NULL) {
    free(napkin->g);
    free(napkin->p);

    stats_free(napkin->tde);

    gsl_rng_free(napkin->rng);
  }

  free(napkin);
}

__attribute__ ((__malloc__))
static struct napkin* napkin_alloc(size_t const ndim,
    size_t const npoly, size_t const nbead,
    size_t const nsubdiv,
    size_t const nthrm, size_t const nprod,
    size_t const nthrmrec, size_t const nprodrec) {
  bool p = true;

  struct napkin* const napkin = malloc(sizeof *napkin);
  if (napkin == NULL)
    p = false;
  else {
    napkin->rng = gsl_rng_alloc(gsl_rng_env_setup());
    if (napkin->rng == NULL)
      p = false;
    else if (!gen_id(napkin->rng, napkin->id, sizeof napkin->id))
      p = false;

    napkin->ensem.Nmemb.dim = ndim;
    napkin->ensem.Nmemb.poly = npoly;
    napkin->ensem.Nmemb.bead = nbead;
    napkin->ensem.Nmemb.subdiv = nsubdiv;

    napkin->ensem.Nstep.thrm = nthrm;
    napkin->ensem.Nstep.prod = nprod;
    napkin->ensem.Nstep.thrmrec = nthrmrec;
    napkin->ensem.Nstep.prodrec = nprodrec;

    napkin->ensem.istep.thrm = 0;
    napkin->ensem.istep.prod = 0;
    napkin->ensem.istep.thrmrec = 0;
    napkin->ensem.istep.prodrec = 0;

    napkin->tde = stats_alloc(nprod);
    if (napkin->tde == NULL)
      p = false;

    size_t const Nbin = size_pow(nsubdiv, ndim);
    napkin->p = malloc(Nbin * sizeof *napkin->p);
    if (napkin->p == NULL)
      p = false;
    else
      for (size_t ibin = 0; ibin < Nbin; ++ibin)
        napkin->p[ibin] = 0;

    napkin->g = malloc(Nbin * sizeof *napkin->g);
    if (napkin->g == NULL)
      p = false;
    else
      for (size_t ibin = 0; ibin < Nbin; ++ibin)
        napkin->g[ibin] = 0;

    napkin->params.ssm.accepted = 0;
    napkin->params.ssm.rejected = 0;
    napkin->params.ssm.dx = 1.0;

    napkin->params.cmd.accepted = 0;
    napkin->params.cmd.rejected = 0;
    napkin->params.cmd.dx = 1.0;
  }

  if (p)
    return napkin;
  else {
    napkin_free(napkin);

    return NULL;
  }
}

void sim_run(size_t const ndim,
    size_t const npoly, size_t const nbead,
    size_t const nsubdiv,
    size_t const nthrm, size_t const nprod,
    size_t const nthrmrec, size_t const nprodrec,
    bool periodic, double L, double beta,
    double (* Vint)(struct ensem const*, struct bead const*, struct bead const*),
    double (* Vend)(struct ensem const*, struct bead const*, struct bead const*),
    double (* Vext)(struct ensem const*, struct bead const*)) {
  err_reset();

  // These assertions guarantee that adding or multiplying two steps or indices
  // never wraps around, which makes their manipulation easier.
  dynamic_assert(ndim <= DIM_MAX, "too many dimensions");
  dynamic_assert(npoly <= POLY_MAX, "too many polymers");
  dynamic_assert(nbead <= BEAD_MAX, "too many beads");
  dynamic_assert(nthrm <= SQRT_SIZE_MAX, "too many thermalization steps");
  dynamic_assert(nprod <= SQRT_SIZE_MAX, "too many production steps");
  dynamic_assert(nthrmrec <= nthrm, "too many thermalization recording steps");
  dynamic_assert(nprodrec <= nprod, "too many production recording steps");

  struct napkin* const napkin = napkin_alloc(ndim, npoly, nbead, nsubdiv,
      nthrm, nprod, nthrmrec, nprodrec);
  if (napkin == NULL)
    err_abort(napkin_alloc);

#ifdef DEBUG

  napkin->ensem.Nstep.thrm = 1;
  napkin->ensem.Nstep.prod = 1;
  napkin->ensem.Nstep.thrmrec = 1;
  napkin->ensem.Nstep.prodrec = 1;

#endif

  napkin->ensem.periodic = periodic;
  napkin->ensem.L = L;
  napkin->ensem.beta = beta;
  napkin->ensem.tau = napkin->ensem.beta / (double) napkin->ensem.Nmemb.bead;
  napkin->ensem.Vint = Vint;
  napkin->ensem.Vend = Vend;
  napkin->ensem.Vext = Vext;

  prepare_stuff(napkin);

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

  Rr_circlatt(napkin);
  Rend_close(napkin);
  // Rend_open(napkin);
  Rm_const(napkin, 1.0);

  save_const(napkin);

  FILE* const paramsfp = res_open(napkin, "params");
  if (paramsfp == NULL)
    err_abort(res_open);

  FILE* const energyfp = res_open(napkin, "energy");
  if (energyfp == NULL)
    err_abort(res_open);

  static double (* const workers[])(struct napkin*) = {work_ss, work_cmd};

  for (size_t istep = 0;
      istep < napkin->ensem.Nstep.thrm + napkin->ensem.Nstep.prod;
      ++istep) {
    choose(napkin, workers
        [ran_index(napkin->rng, sizeof workers / sizeof *workers)]
        (napkin));

    int signum;
    if (sigs_use(&signum))
      switch (signum) {
        case SIGUSR1:
          show_progress(napkin);
          break;
        case SIGUSR2:
          save_mut(napkin);
          break;
      }

    if (istep < napkin->ensem.Nstep.thrm) {
      if (napkin->ensem.Nstep.prodrec * napkin->ensem.istep.thrm >
          napkin->ensem.Nstep.prod * napkin->ensem.istep.thrmrec) {
        disp_params(napkin, paramsfp);

        ++napkin->ensem.istep.thrmrec;
      }

      ++napkin->ensem.istep.thrm;
    } else {
      (void) stats_accum(napkin->tde, est_pimc_tde(&napkin->ensem));

      prob_accum(napkin);
      paircorr_accum(napkin);

      if (napkin->ensem.Nstep.prodrec * napkin->ensem.istep.prod >
          napkin->ensem.Nstep.prod * napkin->ensem.istep.prodrec) {
        disp_params(napkin, paramsfp);
        disp_E(napkin, energyfp);

        ++napkin->ensem.istep.prodrec;
      }

      ++napkin->ensem.istep.prod;
    }
  }

  if (!res_close(napkin, energyfp))
    err_abort(res_close);

  if (!res_close(napkin, paramsfp))
    err_abort(res_close);

  save_mut(napkin);
  show_progress(napkin);
  // show_results(napkin);

  napkin_free(napkin);
}
