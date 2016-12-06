#include "err.h"
#include "exts.h"
#include "fp.h"
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

// These preprocessor directives define
// the limits of the fields in `struct memb`.
#define DIM_MAX ((size_t) 1 << 2) // 4
#define POLY_MAX ((size_t) 1 << 5) // 32
#define BEAD_MAX ((size_t) 1 << 8) // 256
#define SUBDIV_MAX ((size_t) 1 << 11) // 2.048e+3

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

// Vector Math (bead) ---------------------------------------------------------

// Minimal norm squared from origin.
__attribute__ ((__nonnull__, __pure__))
static double bead_norm2(struct ensem const* const ensem,
    struct bead const* const r) {
  double s = 0.0;

  double (* const f)(double, double) =
    ensem->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ensem->Nmemb.dim; ++idim)
    s += gsl_pow_2(f(r->d[idim], ensem->L));

  return s;
}

// Minimal norm from origin.
__attribute__ ((__nonnull__, __pure__))
static double bead_norm(struct ensem const* const ensem,
    struct bead const* const r) {
  return sqrt(bead_norm2(ensem, r));
}

// Minimal distance squared between two points.
__attribute__ ((__nonnull__, __pure__))
static double bead_dist2(struct ensem const* const ensem,
    struct bead const* const r0, struct bead const* const r1) {
  double s = 0.0;

  double (* const f)(double, double) =
    ensem->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ensem->Nmemb.dim; ++idim)
    s += gsl_pow_2(f(r1->d[idim] - r0->d[idim], ensem->L));

  return s;
}

// Minimal distance between two points.
__attribute__ ((__nonnull__, __pure__))
static double bead_dist(struct ensem const* const ensem,
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

// Potential energy for one bead over all polymers.
__attribute__ ((__nonnull__, __pure__))
static double V_bead(struct ensem const* const ensem,
    size_t const ibead) {
  double V = Vint_bead(ensem, ibead) + Vend_bead(ensem, ibead);

  for (size_t ipoly = 0; ipoly < ensem->Nmemb.poly; ++ipoly)
    V += Vext_polybead(ensem, ipoly, ibead);

  return V;
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

// TODO Thermodynamic energy estimator for open polymers.
__attribute__ ((__nonnull__, __pure__))
static double est_pigs_tde(struct ensem const* const ensem) {
  double const E = (double) (ensem->Nmemb.dim * ensem->Nmemb.poly) /
    (2.0 * ensem->tau);
  double const K = K_total(ensem) /
    (double) (ensem->Nmemb.poly * ensem->Nmemb.bead);
  double const V = V_total(ensem) /
    (double) (ensem->Nmemb.poly * ensem->Nmemb.bead);

  return E - K + V;
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
static bool with_res(struct napkin const* const napkin, char const* const str,
    bool (* const f)(struct napkin const*, FILE*)) {
  // Sanity check.
  if (strchr(str, '.') != NULL || strchr(str, '/') != NULL)
    return false;

  char buf[BUFSIZ];
  int const k = snprintf(buf, sizeof buf,
      "%s/%s.data", napkin->id, str);
  if (k < 0 || (size_t) k >= sizeof buf)
    return false;

  FILE* const fp = fopen(buf, "w");
  if (fp == NULL)
    return false;

  bool const p = f(napkin, fp);

  if (fclose(fp) == EOF)
    return false;

  return p;
}

__attribute__ ((__nonnull__))
static bool disp_M(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.Nmemb.bead) < 0)
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
static bool disp_pimc_tde(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu %f %f %f\n",
        napkin->ensem.istep.prod, est_pimc_tde(&napkin->ensem),
        stats_mean(napkin->tde), stats_sem(napkin->tde)) < 0)
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

// Old Crap Printing Stuff (save) ---------------------------------------------

__attribute__ ((__nonnull__))
static void save_Nmemb(struct napkin const* const napkin) {
  char str[BUFSIZ];
  if (snprintf(str, BUFSIZ, "qho-beads-%zud.data", napkin->ensem.Nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  (void) fprintf(fp, "%zu\n", napkin->ensem.Nmemb.bead);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_length(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-length-%zud.data", napkin->ensem.Nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  (void) fprintf(fp, "%f\n", napkin->ensem.L);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_length_harder(struct napkin const* const napkin) {
  char const str[] = "qho-length.data";

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  (void) fprintf(fp, "%f\n", napkin->ensem.L);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_esl(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-ensemble-%zud.data", napkin->ensem.Nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  disp_R(napkin, fp);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_subdivisions(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-subdivisions-%zud.data", napkin->ensem.Nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  (void) fprintf(fp, "%zu\n", napkin->ensem.Nmemb.subdiv);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_potential(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-potential-%zud.data", napkin->ensem.Nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

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
      (void) fprintf(fp, "%f ", r.d[idim]);

    // TODO Unclamp.
    (void) fprintf(fp, "%f %f\n",
        napkin->ensem.Vext(&napkin->ensem, &r),
        fp_clamp(napkin->ensem.Vint(&napkin->ensem, &r0, &r), -8.0, 8.0));
  }

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_probability(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-probability-%zud.data", napkin->ensem.Nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  size_t const Nbin = size_pow(napkin->ensem.Nmemb.subdiv,
      napkin->ensem.Nmemb.dim);

  size_t N = 0.0;
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
      (void) fprintf(fp, "%f ", r.d[idim]);

    (void) fprintf(fp, "%f\n", (double) napkin->p[ibin] / S);
  }

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_paircorr(struct napkin const* const napkin) {
  char const str[] = "qho-paircorrelation.data";

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

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

    (void) fprintf(fp, "%f %f\n", r, (double) napkin->g[ibin] / S);
  }

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

// Other Crap () --------------------------------------------------------------

__attribute__ ((__nonnull__))
static void status_line(struct napkin const* const napkin) {
  (void) printf("i / N = (%zu + %zu) / (%zu + %zu) = %zu %%\n"
      "E = %f +- %f (kappa = 0)\n",
      napkin->ensem.istep.thrm, napkin->ensem.istep.prod, napkin->ensem.Nstep.thrm, napkin->ensem.Nstep.prod,
      100 * (napkin->ensem.istep.thrm + napkin->ensem.istep.prod) /
      (napkin->ensem.Nstep.thrm + napkin->ensem.Nstep.prod),
      stats_mean(napkin->tde), stats_sem(napkin->tde));
}

__attribute__ ((__nonnull__))
static void slow_line(struct napkin const* const napkin) {
  (void) printf("i / N = (%zu + %zu) / (%zu + %zu) = %zu %%\n"
      "E = %f +- %f (kappa = %f)\n",
      napkin->ensem.istep.thrm, napkin->ensem.istep.prod, napkin->ensem.Nstep.thrm, napkin->ensem.Nstep.prod,
      100 * (napkin->ensem.istep.thrm + napkin->ensem.istep.prod) /
      (napkin->ensem.Nstep.thrm + napkin->ensem.Nstep.prod),
      stats_mean(napkin->tde), stats_corrsem(napkin->tde),
      stats_corrtime(napkin->tde));
}

__attribute__ ((__nonnull__))
static void status_file(struct napkin const* const napkin) {
  save_esl(napkin);
  save_probability(napkin);
  save_paircorr(napkin);
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

// More Crap () ---------------------------------------------------------------

__attribute__ ((__const__, __nonnull__))
static double V_zero(
    __attribute__ ((__unused__)) struct ensem const* const ensem,
    __attribute__ ((__unused__)) struct bead const* const r0,
    __attribute__ ((__unused__)) struct bead const* const r1) {
  return 0.0;
}

static double const epsilon = 4.0;
static double const sigma = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double V_lj612(struct ensem const* const ensem,
    struct bead const* const r0, struct bead const* const r1) {
  double const sigmar2 = gsl_pow_2(sigma) / bead_dist2(ensem, r0, r1);

  return 4.0 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

__attribute__ ((__nonnull__, __const__))
static double Vext_zero(
    __attribute__ ((__unused__)) struct ensem const* const ensem,
    __attribute__ ((__unused__)) struct bead const* const r) {
  return 0.0;
}

static double const omega = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double Vext_harm(struct ensem const* const ensem,
    struct bead const* const r) {
  return gsl_pow_2(omega) * bead_norm2(ensem, r) / 2.0;
}

static void simulate(size_t const ndim,
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

  (void) printf("napkin->id = \"%s\"\n", napkin->id);
  prepare_stuff(napkin);

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

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

  Rr_circlatt(napkin);
  Rend_close(napkin);
  // Rend_open(napkin);
  Rm_const(napkin, 1.0);

  /*
  (void) with_res(napkin, "L", disp_L);
  (void) with_res(napkin, "M", disp_M);
  (void) with_res(napkin, "R", disp_R);
  */

  save_length(napkin);
  save_length_harder(napkin);
  save_Nmemb(napkin);
  save_subdivisions(napkin);
  save_potential(napkin);

  FILE* const driftfp = fopen("qho-drift.data", "w");
  if (driftfp == NULL)
    err_abort(fopen);

  FILE* const tdefp = fopen("qho-tde.data", "w");
  if (tdefp == NULL)
    err_abort(fopen);

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
          status_line(napkin);
          break;
        case SIGUSR2:
          status_file(napkin);
          break;
      }

    if (istep < napkin->ensem.Nstep.thrm) {
      if (napkin->ensem.Nstep.prodrec * napkin->ensem.istep.thrm >
          napkin->ensem.Nstep.prod * napkin->ensem.istep.thrmrec) {
        disp_params(napkin, driftfp);

        ++napkin->ensem.istep.thrmrec;
      }

      ++napkin->ensem.istep.thrm;
    } else {
      (void) stats_accum(napkin->tde, est_pimc_tde(&napkin->ensem));

      prob_accum(napkin);
      paircorr_accum(napkin);

      if (napkin->ensem.Nstep.prodrec * napkin->ensem.istep.prod >
          napkin->ensem.Nstep.prod * napkin->ensem.istep.prodrec) {
        disp_params(napkin, driftfp);
        disp_pimc_tde(napkin, tdefp);

        ++napkin->ensem.istep.prodrec;
      }

      ++napkin->ensem.istep.prod;
    }
  }

  if (fclose(tdefp) == EOF)
    err_abort(fclose);

  if (fclose(driftfp) == EOF)
    err_abort(fclose);

  // slow_line(napkin);
  status_line(napkin);
  status_file(napkin);

  napkin_free(napkin);
}

// Main and Procrastination () ------------------------------------------------

// TODO Try other Marsaglia's rngs (UNI or VNI had no observable effect, done).
// TODO Split data files by dimension (use `snprintf`, done).
// TODO Abstract mean/var (done).
// TODO Meditate about the anisotropy and symmetry of potentials (done).
// TODO Deal with the disparity `/ 2` in K_total and V_total (done).
// TODO Make periodicity conditional (done).
// TODO Fix V/K scaling and offsets (done).
// TODO Find home for lost souls (done).
// TODO Consider different masses for different polymers (done).
// TODO Correlation time estimation (done).
// TODO Histograms (position, pair correlation) for the people (done).
// TODO Abstract generic calculations for qho and He-4 (half done).
// TODO PIGS time (half done)!

int main(void) {
  size_t const ndim = 1;
  size_t const npoly = 1;
  size_t const nbead = 32;
  size_t const nsubdiv = ndim == 1 ? 256 : ndim == 2 ? 16 : 4;
  size_t const nthrm = 1 << 14;
  size_t const nprod = 1 << 18;
  size_t const nthrmrec = 1 << 4;
  size_t const nprodrec = 1 << 8;

  double const T = 0.1;
  double const beta = 1.0 / T;

  double const q = omega / 2.0;
  double const E = (double) ndim * q / tanh(q * beta);
  (void) printf("Expected for QHO: E = %f (T = %f)\n", E, T);

  simulate(ndim, npoly, nbead, nsubdiv, nthrm, nprod, nthrmrec, nprodrec,
      true, 10.0, beta, V_lj612, V_zero, Vext_harm);

  return EXIT_SUCCESS;
}
