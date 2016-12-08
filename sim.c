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
      size_t acc;
      size_t rej;
      double h;
      size_t stab;
    } ssm;
    struct {
      size_t acc;
      size_t rej;
      double h;
      size_t stab;
    } cmd;
  } params;
  struct ensem ensem;
  char id[BUFSIZ];
};

double pot_zero(
    __attribute__ ((__unused__)) struct ensem const* const ensem,
    __attribute__ ((__unused__)) struct bead const* const r0,
    __attribute__ ((__unused__)) struct bead const* const r1) {
  return 0.0;
}

double potext_zero(
    __attribute__ ((__unused__)) struct ensem const* const ensem,
    __attribute__ ((__unused__)) struct bead const* const r) {
  return 0.0;
}

double bead_norm2(struct ensem const* const ensem,
    struct bead const* const r) {
  double s = 0.0;

  double (* const f)(double, double) =
    ensem->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ensem->nmemb.dim; ++idim)
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

  for (size_t idim = 0; idim < ensem->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r1->d[idim] - r0->d[idim], ensem->L));

  return s;
}

double bead_dist(struct ensem const* const ensem,
    struct bead const* const r0, struct bead const* const r1) {
  return sqrt(bead_dist2(ensem, r0, r1));
}

// Kinetic energy within a polymer from one bead to the previous one.
static double K_polybead_bw(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  double const M = ensem->R[ipoly].m / (2.0 * gsl_pow_2(ensem->tau));

  if (ibead == 0) {
    size_t const jpoly = ensem->R[ipoly].ifrom;
    size_t const jbead = ensem->nmemb.bead - 1;

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
static double K_polybead_fw(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  double const M = ensem->R[ipoly].m / (2.0 * gsl_pow_2(ensem->tau));

  if (ibead == ensem->nmemb.bead - 1) {
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
static double K_polybead(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  return K_polybead_bw(ensem, ipoly, ibead) +
    K_polybead_fw(ensem, ipoly, ibead);
}

// Kinetic energy within a polymer.
static double K_poly(struct ensem const* const ensem,
    size_t const ipoly) {
  double K = 0.0;

  for (size_t ibead = 0; ibead < ensem->nmemb.bead; ++ibead)
    K += K_polybead_fw(ensem, ipoly, ibead);

  return K;
}

// Kinetic energy.
static double K_total(struct ensem const* const ensem) {
  double K = 0.0;

  for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
    K += K_poly(ensem, ipoly);

  return K;
}

// Potential energy between certain beads of all polymers.
static double Vint_bead(struct ensem const* const ensem,
    size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < ensem->nmemb.poly; ++jpoly)
      V += ensem->Vint(ensem,
          &ensem->R[ipoly].r[ibead], &ensem->R[jpoly].r[ibead]);

  return V;
}

// Additional potential energy between either of the ends of all polymers.
static double Vend_bead(struct ensem const* const ensem,
    size_t const ibead) {
  double V = 0.0;

  if (ibead == 0) {
    for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
      if (ensem->R[ipoly].ifrom == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < ensem->nmemb.poly; ++jpoly)
          if (ensem->R[jpoly].ifrom == SIZE_MAX)
            V += ensem->Vend(ensem,
                &ensem->R[ipoly].r[ibead], &ensem->R[jpoly].r[ibead]);
  } else if (ibead == ensem->nmemb.bead - 1)
    for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
      if (ensem->R[ipoly].ito == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < ensem->nmemb.poly; ++jpoly)
          if (ensem->R[jpoly].ito == SIZE_MAX)
            V += ensem->Vend(ensem,
                &ensem->R[ipoly].r[ibead], &ensem->R[jpoly].r[ibead]);

  return V;
}

// External potential energy for one bead in a polymer.
static double Vext_polybead(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  return ensem->Vext(ensem, &ensem->R[ipoly].r[ibead]);
}

// External potential energy for every bead in a polymer.
static double Vext_bead(struct ensem const* const ensem, size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
    V += Vext_polybead(ensem, ipoly, ibead);

  return V;
}

// Potential energy for one bead over all polymers.
static double V_bead(struct ensem const* const ensem,
    size_t const ibead) {
  return Vint_bead(ensem, ibead) + Vend_bead(ensem, ibead) +
    Vext_bead(ensem, ibead);
}

// Potential energy.
static double V_total(struct ensem const* const ensem) {
  double V = 0.0;

  for (size_t ibead = 0; ibead < ensem->nmemb.bead; ++ibead)
    V += V_bead(ensem, ibead);

  return V;
}

// Thermodynamic energy estimator for closed polymers.
// \langle E_T\rangle & = \frac 1{M \tau} \sum_{k = 1}^M
// \Bigl\langle\frac{d N} 2 - \frac{|R_{k + 1 \bmod M} - R_k|^2}{4 \lambda \tau} + \tau V(R_k)\Bigr\rangle
static double est_pimc_tde(struct ensem const* const ensem) {
  double const E = (double) (ensem->nmemb.dim * ensem->nmemb.poly) /
    (2.0 * ensem->tau);
  double const K = K_total(ensem) /
    (double) (ensem->nmemb.poly * ensem->nmemb.bead);
  double const V = V_total(ensem) /
    (double) (ensem->nmemb.poly * ensem->nmemb.bead);

  return E - K + V;
}

// TODO The crap thermodynamic energy estimator for open polymers.
static double est_pigs_crap(struct ensem const* const ensem) {
  double const V = Vext_bead(ensem, ensem->nmemb.bead / 2) /
    (double) ensem->nmemb.poly;

  return NAN;
}

// TODO Thermodynamic energy estimator for open polymers.
static double est_pigs_tde(struct ensem const* const ensem) {
  return NAN;
}

// Equal mass for every particle.
static void Rm_const(struct napkin* const napkin, double const m) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    napkin->ensem.R[ipoly].m = m;
}

// Close the current configuration.
static void Rend_close(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].ifrom = ipoly;
    napkin->ensem.R[ipoly].ito = ipoly;
  }
}

// Open the current configuration.
static void Rend_open(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].ifrom = SIZE_MAX;
    napkin->ensem.R[ipoly].ito = SIZE_MAX;
  }
}

// Cycle the current configuration.
static void Rend_cycle(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].ifrom =
      size_uwrap_dec(ipoly, napkin->ensem.nmemb.poly);
    napkin->ensem.R[ipoly].ito =
      size_uwrap_inc(ipoly, napkin->ensem.nmemb.poly);
  }
}

// Point initial configuration.
static void Rr_pt(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = fp_lerp(0.5,
            0.0, 1.0, 0.0, napkin->ensem.L);
  }
}

// Random initial configuration.
static void Rr_rand(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = ran_uopen(napkin->rng,
            napkin->ensem.L);
}

// Random point initial configuration.
static void Rr_randpt(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
      double const x = ran_uopen(napkin->rng, napkin->ensem.L);

      for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = x;
    }
}

// Random lattice initial configuration.
static void Rr_randlatt(struct napkin* const napkin) {
  size_t const nlin = size_cirt(napkin->ensem.nmemb.poly,
      napkin->ensem.nmemb.dim);

  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, nlin, i, napkin->ensem.nmemb.dim);

    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = fp_lerp((double) i[idim] +
            ran_uopen(napkin->rng, 1.0),
            0.0, (double) nlin, 0.0, napkin->ensem.L);
  }
}

// Point lattice initial configuration.
static void Rr_ptlatt(struct napkin* const napkin) {
  size_t const nlin = size_cirt(napkin->ensem.nmemb.poly,
      napkin->ensem.nmemb.dim);

  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, nlin, i, napkin->ensem.nmemb.dim);

    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = fp_lerp((double) i[idim],
            0.0, (double) nlin, 0.0, napkin->ensem.L);
  }
}

// Circle lattice initial configuration.
// The equations are based on the theory of Lissajous knots.
static void Rr_circlatt(struct napkin* const napkin) {
  size_t const nlin = size_cirt(napkin->ensem.nmemb.poly,
      napkin->ensem.nmemb.dim);

  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    size_t i[DIM_MAX];
    size_hc(ipoly, nlin, i, napkin->ensem.nmemb.dim);

    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
      double const t = (double) ibead / (double) napkin->ensem.nmemb.bead;

      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
        double const phi = (double) idim / (double) napkin->ensem.nmemb.dim;
        double const x = sin(M_2PI * (t + phi / 2.0));

        napkin->ensem.R[ipoly].r[ibead].d[idim] = fp_lerp((double) i[idim] +
            0.5 + x / 4.0,
            0.0, (double) nlin, 0.0, napkin->ensem.L);
      }
    }
  }
}

static void move_accept_ss(struct napkin* const napkin) {
  ++napkin->params.ssm.acc;
}

static void move_reject_ss(struct napkin* const napkin) {
  size_t const ipoly = napkin->hist.ssm.ipoly;
  size_t const ibead = napkin->hist.ssm.ibead;

  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
    napkin->ensem.R[ipoly].r[ibead].d[idim] = napkin->hist.ssm.r.d[idim];

  ++napkin->params.ssm.rej;
}

static void move_adjust_ss(struct napkin* const napkin) {
  size_t const k = 256;

  napkin->params.ssm.h *= fp_cbalance((double) (napkin->params.ssm.acc + k),
        (double) (napkin->params.ssm.rej + k));
}

static void move_ss(struct napkin* const napkin,
    size_t const ipoly, size_t const ibead) {
  napkin->accept = move_accept_ss;
  napkin->reject = move_reject_ss;
  napkin->adjust = move_adjust_ss;

  napkin->hist.ssm.ipoly = ipoly;
  napkin->hist.ssm.ibead = ibead;

  double (* const f)(double, double) =
    napkin->ensem.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
    napkin->hist.ssm.r.d[idim] = napkin->ensem.R[ipoly].r[ibead].d[idim];

    napkin->ensem.R[ipoly].r[ibead].d[idim] =
      f(napkin->ensem.R[ipoly].r[ibead].d[idim] +
        napkin->params.ssm.h * ran_open(napkin->rng, 1.0), napkin->ensem.L);
  }
}

static void move_accept_cmd(struct napkin* const napkin) {
  ++napkin->params.cmd.acc;
}

static void move_reject_cmd(struct napkin* const napkin) {
  size_t const ipoly = napkin->hist.cmd.ipoly;

  for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
      napkin->ensem.R[ipoly].r[ibead].d[idim] =
        napkin->hist.cmd.R.r[ibead].d[idim];

  ++napkin->params.cmd.rej;
}

static void move_adjust_cmd(struct napkin* const napkin) {
  size_t const k = 256;

  napkin->params.cmd.h *= fp_cbalance((double) (napkin->params.cmd.acc + k),
        (double) (napkin->params.cmd.rej + k));
}

static void move_cmd(struct napkin* const napkin,
    size_t const ipoly) {
  napkin->accept = move_accept_cmd;
  napkin->reject = move_reject_cmd;
  napkin->adjust = move_adjust_cmd;

  napkin->hist.cmd.ipoly = ipoly;

  double (* const f)(double, double) =
    napkin->ensem.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
    double const x = napkin->params.cmd.h *
      ran_open(napkin->rng, 1.0);

    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
      napkin->hist.cmd.R.r[ibead].d[idim] =
        napkin->ensem.R[ipoly].r[ibead].d[idim];

      napkin->ensem.R[ipoly].r[ibead].d[idim] =
        f(napkin->ensem.R[ipoly].r[ibead].d[idim] + x, napkin->ensem.L);
    }
  }
}

static void move_accept_bisect(
    __attribute__ ((__unused__)) struct napkin* const napkin) {
  err_abort(NULL);
}

static void move_reject_bisect(
    __attribute__ ((__unused__)) struct napkin* const napkin) {
  err_abort(NULL);
}

static void move_bisect(
    __attribute__ ((__unused__)) struct napkin* const napkin,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

static void move_accept_swap(
    __attribute__ ((__unused__)) struct napkin* const napkin) {
  err_abort(NULL);
}

static void move_reject_swap(
    __attribute__ ((__unused__)) struct napkin* const napkin) {
  err_abort(NULL);
}

static void move_swap(
    __attribute__ ((__unused__)) struct napkin* const napkin,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

static double work_ss(struct napkin* const napkin) {
  size_t const ipoly = ran_index(napkin->rng, napkin->ensem.nmemb.poly);
  size_t const ibead = ran_index(napkin->rng, napkin->ensem.nmemb.bead);

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

static double work_cmd(struct napkin* const napkin) {
  size_t const ipoly = ran_index(napkin->rng, napkin->ensem.nmemb.poly);

  double const V0 = V_total(&napkin->ensem);

  move_cmd(napkin, ipoly);

  double const V1 = V_total(&napkin->ensem);

  return napkin->ensem.tau * (V1 - V0);
}

static void choose(struct napkin* const napkin,
    double const DeltaS) {
  if (DeltaS > 0.0 && exp(-DeltaS) < ran_uopen(napkin->rng, 1.0))
    napkin->reject(napkin);
  else
    napkin->accept(napkin);

  napkin->adjust(napkin);
}

static void posdist_accum(struct napkin const* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
      size_t i[DIM_MAX];

      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
        i[idim] = size_uclamp(
            (size_t) (floor(fp_lerp(napkin->ensem.R[ipoly].r[ibead].d[idim],
                  0.0, (double) napkin->ensem.L,
                  0.0, (double) napkin->ensem.nmemb.subdiv))),
            napkin->ensem.nmemb.subdiv);

      ++napkin->p[size_unhc(napkin->ensem.nmemb.subdiv, i,
          napkin->ensem.nmemb.dim)];
    }
}

static void paircorr_accum(struct napkin const* const napkin) {
  size_t const nbin = size_pow(napkin->ensem.nmemb.subdiv,
      napkin->ensem.nmemb.dim);

  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < napkin->ensem.nmemb.poly; ++jpoly)
      for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
        size_t const i = size_uclamp(
            (size_t) (floor(fp_lerp(bead_dist(&napkin->ensem,
                    &napkin->ensem.R[ipoly].r[ibead],
                    &napkin->ensem.R[jpoly].r[ibead]),
                  0.0, (double) napkin->ensem.L, 0.0, (double) nbin))), nbin);

        ++napkin->g[i];
      }
}

static bool res_close(
    __attribute__ ((__unused__)) struct napkin const* const napkin,
    FILE* const fp) {
  if (fclose(fp) == EOF)
    return false;

  return true;
}

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

static bool disp_ndim(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.nmemb.dim) < 0)
    return false;

  return true;
}

static bool disp_npoly(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.nmemb.poly) < 0)
    return false;

  return true;
}

static bool disp_nbead(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.nmemb.bead) < 0)
    return false;

  return true;
}

static bool disp_nsubdiv(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu\n", napkin->ensem.nmemb.subdiv) < 0)
    return false;

  return true;
}

static bool disp_length(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%f\n", napkin->ensem.L) < 0)
    return false;

  return true;
}

static bool disp_R_r(struct napkin const* const napkin, FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  if (fprintf(fp, "%zu", iindex) < 0)
    return false;

  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
    if (fprintf(fp, " %f", napkin->ensem.R[ipoly].r[ibead].d[idim]) < 0)
      return false;

  if (fprintf(fp, "\n") < 0)
    return false;

  return true;
}

static bool disp_pots(struct napkin const* const napkin, FILE* const fp) {
  size_t const npt = size_pow(napkin->ensem.nmemb.subdiv + 1,
      napkin->ensem.nmemb.dim);

  struct bead r0;
  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
    r0.d[idim] = 0.0;

  for (size_t ipt = 0; ipt < npt; ++ipt) {
    size_t i[DIM_MAX];
    size_hc(ipt, napkin->ensem.nmemb.subdiv + 1, i, napkin->ensem.nmemb.dim);

    struct bead r;

    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim],
          0.0, (double) napkin->ensem.nmemb.subdiv, 0.0, napkin->ensem.L);

    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
      if (fprintf(fp, "%f ", r.d[idim]) < 0)
        return false;

    if (fprintf(fp, "%f %f\n",
          napkin->ensem.Vext(&napkin->ensem, &r),
          napkin->ensem.Vint(&napkin->ensem, &r0, &r)) < 0)
      return false;
  }

  return true;
}

static bool disp_energy(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu %f %f %f\n",
        napkin->ensem.istep.prod, est_pimc_tde(&napkin->ensem),
        stats_mean(napkin->tde), stats_sem(napkin->tde)) < 0)
    return false;

  return true;
}

static bool disp_params(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "%zu %zu %zu %zu %f %zu %zu %f\n",
        napkin->ensem.istep.thrm, napkin->ensem.istep.prod,
        napkin->params.ssm.acc, napkin->params.ssm.rej,
        napkin->params.ssm.h,
        napkin->params.cmd.acc, napkin->params.cmd.rej,
        napkin->params.cmd.h) < 0)
    return false;

  return true;
}

static bool disp_polys(struct napkin const* const napkin, FILE* const fp) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
      if (!disp_R_r(napkin, fp, ibead, ipoly, ibead))
        return false;

    if (napkin->ensem.R[ipoly].ito != SIZE_MAX)
      if (!disp_R_r(napkin, fp,
            napkin->ensem.nmemb.bead, napkin->ensem.R[ipoly].ito, 0))
        return false;

    if (fprintf(fp, "\n") < 0)
      return false;
  }

  return true;
}

static bool disp_posdist(struct napkin const* const napkin, FILE* const fp) {
  size_t const nbin = size_pow(napkin->ensem.nmemb.subdiv,
      napkin->ensem.nmemb.dim);

  size_t N = 0;
  for (size_t ibin = 0; ibin < nbin; ++ibin)
    N += napkin->p[ibin];
  double const S = (double) N * napkin->ensem.L;

  for (size_t ibin = 0; ibin < nbin; ++ibin) {
    size_t i[DIM_MAX];
    size_hc(ibin, napkin->ensem.nmemb.subdiv, i, napkin->ensem.nmemb.dim);

    struct bead r;

    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
      r.d[idim] = fp_lerp((double) i[idim] + 0.5,
          0.0, (double) napkin->ensem.nmemb.subdiv, 0.0, napkin->ensem.L);

    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
      if (fprintf(fp, "%f ", r.d[idim]) < 0)
        return false;

    if (fprintf(fp, "%f\n", (double) napkin->p[ibin] / S) < 0)
      return false;
  }

  return true;
}

static bool disp_paircorr(struct napkin const* const napkin, FILE* const fp) {
  size_t const nbin = size_pow(napkin->ensem.nmemb.subdiv,
      napkin->ensem.nmemb.dim);

  size_t N = 0.0;
  for (size_t ibin = 0; ibin < nbin; ++ibin)
    N += napkin->g[ibin];
  double const S = N == 0 ? 1.0 : (double) N;

  for (size_t ibin = 0; ibin < nbin; ++ibin) {
    double const r = fp_lerp((double) ibin + 0.5,
        0.0, (double) nbin, 0.0, napkin->ensem.L);

    if (fprintf(fp, "%f %f\n", r, (double) napkin->g[ibin] / S) < 0)
      return false;
  }

  return true;
}

static bool disp_progress(struct napkin const* const napkin, FILE* const fp) {
  size_t const i = napkin->ensem.istep.thrm + napkin->ensem.istep.prod;
  size_t const n = napkin->ensem.nstep.thrm + napkin->ensem.nstep.prod;

  if (fprintf(fp, "i / n = (%zu + %zu) / (%zu + %zu) = %zu %%\n",
      napkin->ensem.istep.thrm, napkin->ensem.istep.prod,
      napkin->ensem.nstep.thrm, napkin->ensem.nstep.prod,
      n == 0 ? 100 : 100 * i / n) < 0)
    return false;

  return true;
}

static bool disp_results(struct napkin const* const napkin, FILE* const fp) {
  if (fprintf(fp, "E = %f +- %f (kappa = %f)\n",
      stats_mean(napkin->tde), stats_corrsem(napkin->tde),
      stats_corrtime(napkin->tde)) < 0)
    return false;

  return true;
}

static bool disp_wrong_results_fast(struct napkin const* const napkin,
    FILE* const fp) {
  if (fprintf(fp, "E = %f +- %f (kappa = %f)\n",
        stats_mean(napkin->tde), 1000.0 * stats_sem(napkin->tde), 1000.0) < 0)
    return false;

  return true;
}

static bool prepare_save(struct napkin const* const napkin) {
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

static bool save_const(struct napkin const* const napkin) {
  return res_use(napkin, "ndim", disp_ndim) &&
    res_use(napkin, "npoly", disp_npoly) &&
    res_use(napkin, "nbead", disp_nbead) &&
    res_use(napkin, "nsubdiv", disp_nsubdiv) &&
    res_use(napkin, "length", disp_length) &&
    res_use(napkin, "pots", disp_pots);
}

static bool save_mut(struct napkin const* const napkin) {
  return res_use(napkin, "polys", disp_polys) &&
    res_use(napkin, "posdist", disp_posdist) &&
    res_use(napkin, "paircorr", disp_paircorr);
}

static void napkin_free(struct napkin* const napkin) {
  if (napkin != NULL) {
    free(napkin->g);
    free(napkin->p);

    stats_free(napkin->tde);

    gsl_rng_free(napkin->rng);
  }

  free(napkin);
}

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
    else if (!ran_dateid(napkin->rng, "run", napkin->id, sizeof napkin->id))
      p = false;

    napkin->ensem.nmemb.dim = ndim;
    napkin->ensem.nmemb.poly = npoly;
    napkin->ensem.nmemb.bead = nbead;
    napkin->ensem.nmemb.subdiv = nsubdiv;

    napkin->ensem.nstep.thrm = nthrm;
    napkin->ensem.nstep.prod = nprod;
    napkin->ensem.nstep.thrmrec = nthrmrec;
    napkin->ensem.nstep.prodrec = nprodrec;

    napkin->ensem.istep.thrm = 0;
    napkin->ensem.istep.prod = 0;
    napkin->ensem.istep.thrmrec = 0;
    napkin->ensem.istep.prodrec = 0;

    napkin->tde = stats_alloc(nprod);
    if (napkin->tde == NULL)
      p = false;

    size_t const nbin = size_pow(nsubdiv, ndim);
    napkin->p = malloc(nbin * sizeof *napkin->p);
    if (napkin->p == NULL)
      p = false;
    else
      for (size_t ibin = 0; ibin < nbin; ++ibin)
        napkin->p[ibin] = 0.0;

    napkin->g = malloc(nbin * sizeof *napkin->g);
    if (napkin->g == NULL)
      p = false;
    else
      for (size_t ibin = 0; ibin < nbin; ++ibin)
        napkin->g[ibin] = 0.0;

    // TODO Use `napkin->L / 2`.
    napkin->params.ssm.acc = 0;
    napkin->params.ssm.rej = 0;
    napkin->params.ssm.h = 1.0;

    // TODO Use `napkin->L / 2`.
    napkin->params.cmd.acc = 0;
    napkin->params.cmd.rej = 0;
    napkin->params.cmd.h = 1.0;
  }

  if (p)
    return napkin;
  else {
    napkin_free(napkin);

    return NULL;
  }
}

// TODO Should probably split napkin allocation and initialization and
// expose related procedures for doing just that.
bool sim_run(size_t const ndim, size_t const npoly,
    size_t const nbead, size_t const nsubdiv,
    size_t const nthrm, size_t const nprod,
    size_t const nthrmrec, size_t const nprodrec,
    bool const periodic, double const L, double const beta,
    double (* const Vint)(struct ensem const*,
      struct bead const*, struct bead const*),
    double (* const Vend)(struct ensem const*,
      struct bead const*, struct bead const*),
    double (* const Vext)(struct ensem const*, struct bead const*)) {
  err_reset();

#ifdef _GNU_SOURCE
#ifdef DEBUG
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
#endif

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

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

  napkin->ensem.periodic = periodic;
  napkin->ensem.L = L;
  napkin->ensem.beta = beta;
  napkin->ensem.tau = napkin->ensem.beta / (double) napkin->ensem.nmemb.bead;
  napkin->ensem.Vint = Vint;
  napkin->ensem.Vend = Vend;
  napkin->ensem.Vext = Vext;

  Rr_pt(napkin);
  Rend_close(napkin);
  // Rend_open(napkin);
  Rm_const(napkin, 1.0);

  if (!prepare_save(napkin))
    err_abort(prepare_save);

  if (!save_const(napkin))
    err_abort(save_const);

  FILE* const paramsfp = res_open(napkin, "params");
  if (paramsfp == NULL)
    err_abort(res_open);

  FILE* const energyfp = res_open(napkin, "energy");
  if (energyfp == NULL)
    err_abort(res_open);

  static double (* const workers[])(struct napkin*) = {work_ss, work_cmd};

  for (size_t istep = 0;
      istep < napkin->ensem.nstep.thrm + napkin->ensem.nstep.prod;
      ++istep) {
    choose(napkin, workers
        [ran_index(napkin->rng, sizeof workers / sizeof *workers)]
        (napkin));

    int signum;
    if (sigs_use(&signum))
      switch (signum) {
        case SIGUSR1:
          (void) disp_progress(napkin, stdout);
          break;
        case SIGUSR2:
          (void) save_mut(napkin);
          break;
      }

    if (istep < napkin->ensem.nstep.thrm) {
      if (napkin->ensem.nstep.prodrec * napkin->ensem.istep.thrm >
          napkin->ensem.nstep.prod * napkin->ensem.istep.thrmrec) {
        if (!disp_params(napkin, paramsfp))
          err_abort(disp_params);

        ++napkin->ensem.istep.thrmrec;
      }

      ++napkin->ensem.istep.thrm;
    } else {
      (void) stats_accum(napkin->tde, est_pimc_tde(&napkin->ensem));

      posdist_accum(napkin);
      paircorr_accum(napkin);

      if (napkin->ensem.nstep.prodrec * napkin->ensem.istep.prod >
          napkin->ensem.nstep.prod * napkin->ensem.istep.prodrec) {
        if (!disp_energy(napkin, energyfp))
          err_abort(disp_energy);

        if (!disp_params(napkin, paramsfp))
          err_abort(disp_params);

        ++napkin->ensem.istep.prodrec;
      }

      ++napkin->ensem.istep.prod;
    }
  }

  if (!res_close(napkin, energyfp))
    err_abort(res_close);

  if (!res_close(napkin, paramsfp))
    err_abort(res_close);

  if (!save_mut(napkin))
    err_abort(save_mut);

  if (!disp_progress(napkin, stdout))
    err_abort(disp_progress);

  if (!disp_wrong_results_fast(napkin, stdout))
    err_abort(disp_wrong_results_fast);

  napkin_free(napkin);

  return true;
}
