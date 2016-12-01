#include "err.h"
#include "exts.h"
#include "fp.h"
#include "sigs.h"
#include "size.h"
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Types () -------------------------------------------------------------------

struct par {
  size_t accepted;
  size_t rejected;
  double dx;
};

struct est {
  struct {
    size_t N; // Samples.
    double M1; // First moment.
    double M2; // Second moment.
  } acc; // Accumulator.
  double mean; // Mean.
  double var; // Variance.
  double sd; // Standard deviation.
  double sem; // Standard error of the mean.
};

struct step {
  size_t thrm; // Thermalization steps.
  size_t prod; // Production steps.
  size_t thrmrec; // Recording steps during thermalization.
  size_t prodrec; // Recording steps during production.
};

struct memb {
  size_t poly;
  size_t bead;
  size_t dim;
  size_t subdiv;
};

struct bead {
  double* d;
};

struct poly {
  struct bead* r;
  size_t from;
  size_t to;
  double m;
};

struct napkin {
  struct memb nmemb;
  struct step nstep;
  struct step istep;
  gsl_rng* rng; // State of the random number generator.
  struct poly* R;
  struct est* tde; // Thermodynamic energy estimator.
  bool periodic;
  double L; // Lattice constant.
  double beta; // Inverse temperature sometimes.
  double tau; // Imaginary time step maybe.
  double (* Vint)(struct napkin const*, struct bead, struct bead); // Potential between beads.
  double (* Vend)(struct napkin const*, struct bead, struct bead); // Extra potential between open ends.
  double (* Vext)(struct napkin const*, struct bead); // External potential.
  void (* accept)(struct napkin*);
  void (* reject)(struct napkin*);
  void (* adjust)(struct napkin*);
  struct {
    struct {
      size_t ipoly;
      size_t ibead;
      struct bead r;
    } ssm; // Single slice move.
    struct {
      size_t ipoly;
      struct poly R;
    } comd; // Center of mass displacement.
  } hist; // Undo history for the latest moves.
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
    } comd;
  } params; // Optimization parameters for each move.
};

// Lost Souls () --------------------------------------------------------------

static size_t ran_index(gsl_rng* const rng, size_t const n) {
  return (size_t) gsl_rng_uniform_int(rng, n);
}

// Estimator Statistics (est) -------------------------------------------------

static void est_free(struct est* const acc) {
  free(acc);
}

static struct est* est_alloc(void) {
  struct est* const est = malloc(sizeof *est);
  if (est == NULL)
    return NULL;

  est->acc.N = 0;
  est->acc.M1 = 0.0;
  est->acc.M2 = 0.0;

  return est;
}

static void est_accum(struct est* const est, double const x) {
  ++est->acc.N;
  double const dx = x - est->acc.M1;
  est->acc.M1 += dx / (double) est->acc.N;
  est->acc.M2 += dx * (x - est->acc.M1);
}

static void est_calc(struct est* const est) {
  switch (est->acc.N) {
    case 0:
      est->mean = NAN;
      est->var = NAN;
      est->sd = NAN;
      est->sem = NAN;
      break;
    case 1:
      est->mean = est->acc.M1;
      est->var = 0.0;
      est->sd = 0.0;
      est->sem = 0.0;
      break;
    default:
      est->mean = est->acc.M1;
      est->var = est->acc.M2 / (double) (est->acc.N - 1);
      est->sd = sqrt(est->var);
      est->sem = est->sd / sqrt((double) est->acc.N);
  }
}

// Vector Math (bead) ---------------------------------------------------------

// Minimal norm squared from origin.
static double bead_norm2(struct napkin const* const napkin,
    struct bead const r) {
  double s = 0.0;

  double (* const f)(double, double) = napkin->periodic ?
    fp_wrap : fp_constant;

  for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r.d[idim], napkin->L));

  return s;
}

// Minimal norm from origin.
static double bead_norm(struct napkin const* const napkin,
    struct bead const r) {
  return sqrt(bead_norm2(napkin, r));
}

// Minimal distance squared between two points.
static double bead_dist2(struct napkin const* const napkin,
    struct bead const r0, struct bead const r1) {
  double s = 0.0;

  double (* const f)(double, double) = napkin->periodic ?
    fp_wrap : fp_constant;

  for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r1.d[idim] - r0.d[idim], napkin->L));

  return s;
}

// Minimal distance between two points.
static double bead_dist(struct napkin const* const napkin,
    struct bead const r0, struct bead const r1) {
  return sqrt(bead_dist2(napkin, r0, r1));
}

// Mass Assignment (mass) -----------------------------------------------------

// Equal mass for every particle.
static void mass_const(struct napkin* const napkin, double const m) {
  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    napkin->R[ipoly].m = m;
}

// Configuration Forcing (force) ----------------------------------------------

// Close the current configuration.
static void force_close(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly) {
    napkin->R[ipoly].from = ipoly;
    napkin->R[ipoly].to = ipoly;
  }
}

// Open the current configuration.
static void force_open(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly) {
    napkin->R[ipoly].from = SIZE_MAX;
    napkin->R[ipoly].to = SIZE_MAX;
  }
}

// Cycle the current configuration.
static void force_cycle(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly) {
    napkin->R[ipoly].from = size_uwrap(ipoly + napkin->nmemb.poly - 1, napkin->nmemb.poly);
    napkin->R[ipoly].to = size_uwrap(ipoly + 1, napkin->nmemb.poly);
  }
}

// Initial Configurations (conf) ----------------------------------------------

// Random initial configuration.
static void conf_rand(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim)
        napkin->R[ipoly].r[ibead].d[idim] =
          gsl_ran_flat(napkin->rng, 0.0, napkin->L);
}

// Random point initial configuration.
static void conf_randpt(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim) {
      double const x = gsl_ran_flat(napkin->rng, 0, napkin->L);

      for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead)
        napkin->R[ipoly].r[ibead].d[idim] = x;
    }
}

// Random lattice initial configuration.
static void conf_randlatt(struct napkin* const napkin) {
  double const w = ceil(pow(napkin->nmemb.poly, 1.0 / napkin->nmemb.dim));
  double const v = napkin->L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead) {
      size_t i = ipoly;

      for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim) {
        size_div_t const z = size_div(i, n);

        napkin->R[ipoly].r[ibead].d[idim] =
          (double) z.rem * v + gsl_ran_flat(napkin->rng, 0, v);
        i = z.quot;
      }
    }
}

// Point lattice initial configuration.
static void conf_ptlatt(struct napkin* const napkin) {
  double const w = ceil(pow(napkin->nmemb.poly, 1.0 / napkin->nmemb.dim));
  double const v = napkin->L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead) {
      size_t i = ipoly;

      for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim) {
        size_div_t const z = size_div(i, n);

        napkin->R[ipoly].r[ibead].d[idim] = (double) z.rem * v;
        i = z.quot;
      }
    }
}

// Circle lattice initial configuration.
// The equations are based on the theory of Lissajous knots.
static void conf_circlatt(struct napkin* const napkin) {
  double const w = ceil(pow(napkin->nmemb.poly, 1.0 / napkin->nmemb.dim));
  double const v = napkin->L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead) {
      double const t = (double) ibead / napkin->nmemb.bead;

      size_t i = ipoly;

      for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim) {
        double const phi = (double) idim / napkin->nmemb.dim;
        double const y = sin(M_2PI * (t + phi / 2.0)) / 2.0;

        size_div_t const z = size_div(i, n);
        napkin->R[ipoly].r[ibead].d[idim] =
          ((double) z.rem + (y + 1.0) / 2.0) * v;
        i = z.quot;
      }
    }
}

// Physical Moves (move) ------------------------------------------------------

// r <=> a: z >=< 1
// r = 0: z = 3 / 2
// a = 0: z = 1 / 2
static double surface(double const a, double const r) {
  return 0.5 + a / (a + r);
}

// r <=> a: z >=< 1
// r = 0: z = 1 + exp(-a)
// a = 0: z = 1 - exp(-r)
static double another_surface(double const a, double const r) {
  return 1.0 - exp(-a) + exp(-r);
}

static void move_accept_ss(struct napkin* const napkin) {
  ++napkin->params.ssm.accepted;
}

static void move_reject_ss(struct napkin* const napkin) {
  size_t const ipoly = napkin->hist.ssm.ipoly;
  size_t const ibead = napkin->hist.ssm.ibead;

  for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim)
    napkin->R[ipoly].r[ibead].d[idim] = napkin->hist.ssm.r.d[idim];

  ++napkin->params.ssm.rejected;
}

static void move_adjust_ss(struct napkin* const napkin) {
  napkin->params.ssm.dx = napkin->params.ssm.dx *
    surface((double) napkin->params.ssm.accepted,
        (double) napkin->params.ssm.rejected);
}

static void move_ss(struct napkin* const napkin,
    size_t const ipoly, size_t const ibead) {
  napkin->accept = move_accept_ss;
  napkin->reject = move_reject_ss;
  napkin->adjust = move_adjust_ss;

  napkin->hist.ssm.ipoly = ipoly;
  napkin->hist.ssm.ibead = ibead;

  double (* const f)(double, double) = napkin->periodic ?
    fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim) {
    napkin->hist.ssm.r.d[idim] = napkin->R[ipoly].r[ibead].d[idim];

    napkin->R[ipoly].r[ibead].d[idim] = f(
        napkin->R[ipoly].r[ibead].d[idim] +
        napkin->params.ssm.dx * gsl_ran_flat(napkin->rng, -1.0, 1.0),
        napkin->L);
  }
}

static void move_accept_comd(struct napkin* const napkin) {
  ++napkin->params.comd.accepted;
}

static void move_reject_comd(struct napkin* const napkin) {
  size_t const ipoly = napkin->hist.comd.ipoly;

  for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim)
      napkin->R[ipoly].r[ibead].d[idim] =
        napkin->hist.comd.R.r[ibead].d[idim];

  ++napkin->params.comd.rejected;
}

static void move_adjust_comd(struct napkin* const napkin) {
  napkin->params.comd.dx = napkin->params.comd.dx *
    another_surface((double) napkin->params.comd.accepted,
        (double) napkin->params.comd.rejected);
}

static void move_comd(struct napkin* const napkin,
    size_t const ipoly) {
  napkin->accept = move_accept_comd;
  napkin->reject = move_reject_comd;
  napkin->adjust = move_adjust_comd;

  napkin->hist.comd.ipoly = ipoly;

  double (* const f)(double, double) = napkin->periodic ?
    fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim) {
    double const x = napkin->params.comd.dx *
      gsl_ran_flat(napkin->rng, -1.0, 1.0);

    for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead) {
      napkin->hist.comd.R.r[ibead].d[idim] =
        napkin->R[ipoly].r[ibead].d[idim];

      napkin->R[ipoly].r[ibead].d[idim] = f(
          napkin->R[ipoly].r[ibead].d[idim] + x,
          napkin->L);
    }
  }
}

__attribute__ ((__noreturn__))
static void move_accept_bisect(__attribute__ ((__unused__))
    struct napkin* const napkin) {
  err_abort(NULL);
}

__attribute__ ((__noreturn__))
static void move_reject_bisect(__attribute__ ((__unused__))
    struct napkin* const napkin) {
  err_abort(NULL);
}

__attribute__ ((__noreturn__))
static void move_bisect(
    __attribute__ ((__unused__)) struct napkin* const napkin,
    __attribute__ ((__unused__)) size_t const ipoly,
    __attribute__ ((__unused__)) size_t const ibead) {
  err_abort(NULL);
}

// Kinetic and Potential Energies (K, V) --------------------------------------

static double K_polybead_bw(struct napkin const* const napkin,
    size_t const ipoly, size_t const ibead) {
  double const M = napkin->R[ipoly].m /
    (2.0 * gsl_pow_2(napkin->tau));

  if (ibead == 0) {
    size_t const jpoly = napkin->R[ipoly].from;
    size_t const jbead = napkin->nmemb.bead - 1;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return M * bead_dist2(napkin,
          napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead - 1;

    return M * bead_dist2(napkin,
        napkin->R[ipoly].r[ibead], napkin->R[ipoly].r[jbead]);
  }
}

// \lambda = \hbar / 2 m, K = x^2 / 4 \lambda \tau^2 = M x^2
// M = m / 2 \hbar \tau^2
static double K_polybead_fw(struct napkin const* const napkin,
    size_t const ipoly, size_t const ibead) {
  double const M = napkin->R[ipoly].m /
    (2.0 * gsl_pow_2(napkin->tau));

  if (ibead == napkin->nmemb.bead - 1) {
    size_t const jpoly = napkin->R[ipoly].to;
    size_t const jbead = 0;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return M * bead_dist2(napkin,
          napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead + 1;

    return M * bead_dist2(napkin,
        napkin->R[ipoly].r[ibead], napkin->R[ipoly].r[jbead]);
  }
}

static double K_polybead(struct napkin const* const napkin,
    size_t const ipoly, size_t const ibead) {
  return K_polybead_bw(napkin, ipoly, ibead) +
    K_polybead_fw(napkin, ipoly, ibead);
}

static double K_poly(struct napkin const* const napkin,
    size_t const ipoly) {
  double K = 0.0;

  for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead)
    K += K_polybead_fw(napkin, ipoly, ibead);

  return K;
}

static double K_total(struct napkin const* const napkin) {
  double K = 0.0;

  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    K += K_poly(napkin, ipoly);

  return K;
}

static double Vint_bead(struct napkin const* const napkin,
    size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < napkin->nmemb.poly; ++jpoly)
      V += napkin->Vint(napkin,
          napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[ibead]);

  return V;
}

static double Vend_bead(struct napkin const* const napkin,
    size_t const ibead) {
  double V = 0.0;

  if (ibead == 0) {
    for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
      if (napkin->R[ipoly].from == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < napkin->nmemb.poly; ++jpoly)
          if (napkin->R[jpoly].from == SIZE_MAX)
            V += napkin->Vend(napkin,
                napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[ibead]);
  } else if (ibead == napkin->nmemb.bead - 1)
    for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
      if (napkin->R[ipoly].to == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < napkin->nmemb.poly; ++jpoly)
          if (napkin->R[jpoly].to == SIZE_MAX)
            V += napkin->Vend(napkin,
                napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[ibead]);

  return V;
}

static double Vext_polybead(struct napkin const* const napkin,
    size_t const ipoly, size_t const ibead) {
  return napkin->Vext(napkin, napkin->R[ipoly].r[ibead]);
}

static double V_bead(struct napkin const* const napkin,
    size_t const ibead) {
  double V = Vint_bead(napkin, ibead) + Vend_bead(napkin, ibead);

  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly)
    V += Vext_polybead(napkin, ipoly, ibead);

  return V;
}

static double V_total(struct napkin const* const napkin) {
  double V = 0.0;

  for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead)
    V += V_bead(napkin, ibead);

  return V;
}

// Action Stuff (S) -----------------------------------------------------------

static double SV(struct napkin const* const napkin, double const V) {
  return napkin->tau * V;
}

static double SK(struct napkin const* const napkin, double const K) {
  return napkin->tau * K;
}

// Estimators (est) -----------------------------------------------------------

// \langle E_T\rangle & = \frac 1{M \tau} \sum_{k = 1}^M
// \Bigl\langle\frac{d N} 2 - \frac{|R_{k + 1 \bmod M} - R_k|^2}{4 \lambda \tau} + \tau V(R_k)\Bigr\rangle
static double est_tde(struct napkin const* const napkin) {
  double const E = napkin->nmemb.dim * napkin->nmemb.poly /
    (2.0 * napkin->tau);
  double const K = K_total(napkin) / (napkin->nmemb.poly * napkin->nmemb.bead);
  double const V = V_total(napkin) / (napkin->nmemb.poly * napkin->nmemb.bead);

  return E - K + V;
}

// Workers (work) -------------------------------------------------------------

static double work_ss(struct napkin* const napkin) {
  size_t const ipoly = ran_index(napkin->rng, napkin->nmemb.poly);
  size_t const ibead = ran_index(napkin->rng, napkin->nmemb.bead);

  double const V0 =
    Vint_bead(napkin, ibead) + Vend_bead(napkin, ibead) +
    Vext_polybead(napkin, ipoly, ibead);
  double const K0 = K_polybead(napkin, ipoly, ibead);

  move_ss(napkin, ipoly, ibead);

  double const V1 =
    Vint_bead(napkin, ibead) + Vend_bead(napkin, ibead) +
    Vext_polybead(napkin, ipoly, ibead);
  double const K1 = K_polybead(napkin, ipoly, ibead);

  return SK(napkin, K1 - K0) + SV(napkin, V1 - V0);
}

static double work_comd(struct napkin* const napkin) {
  size_t const ipoly = ran_index(napkin->rng, napkin->nmemb.poly);

  double const V0 = V_total(napkin);

  move_comd(napkin, ipoly);

  double const V1 = V_total(napkin);

  return SV(napkin, V1 - V0);
}

// Work Horses () -------------------------------------------------------------

static void choose(struct napkin* const napkin,
    double const DeltaS) {
  if (DeltaS <= 0.0 || gsl_ran_flat(napkin->rng, 0.0, 1.0) <= exp(-DeltaS))
    napkin->accept(napkin);
  else
    napkin->reject(napkin);

  napkin->adjust(napkin);
}

// Printing into Data Files (disp) --------------------------------------------

static void disp_bead(struct napkin const* const napkin,
    FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  (void) fprintf(fp, "%zu", iindex);

  for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim)
    (void) fprintf(fp, " %f", napkin->R[ipoly].r[ibead].d[idim]);

  (void) fprintf(fp, "\n");
}

static void disp_poly(struct napkin const* const napkin,
    FILE* const fp) {
  for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly) {
    for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead)
      disp_bead(napkin, fp, ibead, ipoly, ibead);

    if (napkin->R[ipoly].to != SIZE_MAX)
      disp_bead(napkin, fp, napkin->nmemb.bead, napkin->R[ipoly].to, 0);

    (void) fprintf(fp, "\n");
  }
}

static void disp_tde(struct napkin const* const napkin,
    FILE* const fp) {
  est_calc(napkin->tde);

  (void) fprintf(fp, "%zu %f %f %f\n",
      // Correlations!
      napkin->istep.prod, est_tde(napkin),
      napkin->tde->mean, napkin->tde->sem * napkin->nmemb.bead);
}

static void disp_drift(struct napkin const* const napkin,
    FILE* const fp) {
  (void) fprintf(fp, "%zu %f %f\n",
      napkin->istep.thrm + napkin->istep.prod,
      napkin->params.ssm.dx, napkin->params.comd.dx);
}

// Saving into Data Files (save)

static void save_length(struct napkin const* const napkin) {
  FILE* const fp = fopen("qho-length.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  (void) fprintf(fp, "%f\n", napkin->L);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

static void save_esl(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-ensemble-%zud.data", napkin->nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  disp_poly(napkin, fp);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

static void save_subdivisions(struct napkin const* const napkin) {
  FILE* const fp = fopen("qho-subdivisions.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  (void) fprintf(fp, "%zu\n", napkin->nmemb.subdiv);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

static void save_potential(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-potential-%zud.data", napkin->nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  double const v = napkin->L / napkin->nmemb.subdiv;

  double* const d = malloc(napkin->nmemb.dim * sizeof *d);
  if (d == NULL)
    err_abort(malloc);

  struct bead r = {.d = d};

  for (size_t ipt = 0; ipt < size_pow(napkin->nmemb.subdiv + 1, napkin->nmemb.dim); ++ipt) {
    size_t i = ipt;

    for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim) {
      size_div_t const z = size_div(i, napkin->nmemb.subdiv + 1);

      r.d[idim] = (double) z.rem * v;
      i = z.quot;
    }

    for (size_t idim = 0; idim < napkin->nmemb.dim; ++idim)
      (void) fprintf(fp, "%f ", r.d[idim]);

    (void) fprintf(fp, "%f\n", napkin->Vext(napkin, r));
  }

  free(d);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

// Other Crap () --------------------------------------------------------------

static void status_line(struct napkin const* const napkin) {
  est_calc(napkin->tde);

  (void) printf("i / N = (%zu + %zu) / (%zu + %zu) = %zu %%\n"
      "E = %f +- %f\n",
      napkin->istep.thrm, napkin->istep.prod, napkin->nstep.thrm, napkin->nstep.prod,
      100 *
      (napkin->istep.thrm + napkin->istep.prod) / (napkin->nstep.thrm + napkin->nstep.prod),
      // Correlations!
      napkin->tde->mean, napkin->tde->sem * napkin->nmemb.bead);
}

static void status_file(struct napkin const* const napkin) {
  save_esl(napkin);
}

// More Crap () ---------------------------------------------------------------

static double const epsilon = 1.0;
static double const sigma = 1.0;

static double lj612b(struct napkin const* const napkin,
    struct bead const r0, struct bead const r1) {
  double const sigmar2 = gsl_pow_2(sigma) / bead_dist2(napkin, r0, r1);

  return 4.0 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

static double const omega = 1.0;

static double harmb(struct napkin const* const napkin,
    struct bead const r) {
  return napkin->R[0].m * gsl_pow_2(omega) * bead_norm2(napkin, r) / 2.0;
}

static double zerob(__attribute__ ((__unused__))
    struct napkin const* const napkin,
    __attribute__ ((__unused__)) struct bead const r0,
    __attribute__ ((__unused__)) struct bead const r1) {
  return 0.0;
}

// Allocation (alloc) ---------------------------------------------------------

typedef void* (* alloc)(size_t);

static void* alloc_err(__attribute__ ((__unused__)) size_t const n) {
  return NULL;
}

// Napkin Management (napkin) -------------------------------------------------

static void napkin_free(struct napkin* const napkin) {
  if (napkin != NULL) {
    if (napkin->hist.comd.R.r != NULL)
      for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead)
        free(napkin->hist.comd.R.r[ibead].d);

    free(napkin->hist.comd.R.r);

    free(napkin->hist.ssm.r.d);

    est_free(napkin->tde);

    if (napkin->R != NULL)
      for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly) {
        if (napkin->R[ipoly].r != NULL)
          for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead)
            free(napkin->R[ipoly].r[ibead].d);

        free(napkin->R[ipoly].r);
      }

    free(napkin->R);

    gsl_rng_free(napkin->rng);
  }

  free(napkin);
}

static struct napkin* napkin_alloc(size_t const npoly, size_t const nbead,
    size_t const ndim, size_t const nsubdiv,
    size_t const nthrm, size_t const nprod,
    size_t const nthrmrec, size_t const nprodrec) {
  // These assertions guarantee that adding or multiplying two steps or indices
  // never wraps around, which makes their manipulation easier.
  dynamic_assert(npoly <= RT_SIZE_MAX(4), "too many polymers");
  dynamic_assert(nbead <= RT_SIZE_MAX(4), "too many beads");
  dynamic_assert(ndim <= RT_SIZE_MAX(4), "too many dimensions");
  dynamic_assert(nthrm <= RT_SIZE_MAX(2), "too many thermalization steps");
  dynamic_assert(nprod <= RT_SIZE_MAX(2), "too many production steps");
  dynamic_assert(nthrmrec <= nthrm, "too many thermalization recording steps");
  dynamic_assert(nprodrec <= nprod, "too many production recording steps");

  alloc alloc = malloc;

  struct napkin* const napkin = alloc(sizeof *napkin);
  if (napkin == NULL)
    alloc = alloc_err;
  else {
    napkin->nmemb.poly = npoly;
    napkin->nmemb.bead = nbead;
    napkin->nmemb.dim = ndim;
    napkin->nmemb.subdiv = nsubdiv;

    napkin->nstep.thrm = nthrm;
    napkin->nstep.prod = nprod;
    napkin->nstep.thrmrec = nthrm;
    napkin->nstep.prodrec = nprodrec;

    napkin->istep.thrm = 0;
    napkin->istep.prod = 0;
    napkin->istep.thrmrec = 0;
    napkin->istep.prodrec = 0;

    gsl_rng_type const* const t = gsl_rng_env_setup();
    napkin->rng = gsl_rng_alloc(t);
    if (napkin->rng == NULL)
      alloc = alloc_err;

    napkin->R = alloc(napkin->nmemb.poly * sizeof *napkin->R);
    if (napkin->R == NULL)
      alloc = alloc_err;
    else
      for (size_t ipoly = 0; ipoly < napkin->nmemb.poly; ++ipoly) {
        napkin->R[ipoly].r =
          alloc(napkin->nmemb.bead * sizeof *napkin->R[ipoly].r);
        if (napkin->R[ipoly].r == NULL)
          alloc = alloc_err;
        else
          for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead) {
            napkin->R[ipoly].r[ibead].d =
              alloc(napkin->nmemb.dim * sizeof *napkin->R[ipoly].r[ibead].d);
            if (napkin->R[ipoly].r[ibead].d == NULL)
              alloc = alloc_err;
          }
      }

    napkin->hist.ssm.r.d =
      alloc(napkin->nmemb.dim * sizeof *napkin->hist.ssm.r.d);
    if (napkin->hist.ssm.r.d == NULL)
      alloc = alloc_err;

    napkin->hist.comd.R.r =
      alloc(napkin->nmemb.bead * sizeof *napkin->hist.comd.R.r);
    if (napkin->hist.comd.R.r == NULL)
      alloc = alloc_err;
    else
      for (size_t ibead = 0; ibead < napkin->nmemb.bead; ++ibead) {
        napkin->hist.comd.R.r[ibead].d =
          alloc(napkin->nmemb.dim * sizeof *napkin->hist.comd.R.r[ibead].d);
        if (napkin->hist.comd.R.r[ibead].d == NULL)
          alloc = alloc_err;
      }

    napkin->tde = est_alloc();
    if (napkin->tde == NULL)
      alloc = alloc_err;

    napkin->params.ssm.accepted = 0;
    napkin->params.ssm.rejected = 0;
    napkin->params.ssm.dx = 1;

    napkin->params.comd.accepted = 0;
    napkin->params.comd.rejected = 0;
    napkin->params.comd.dx = 1;
  }

  if (alloc == alloc_err) {
    napkin_free(napkin);

    return NULL;
  } else
    return napkin;
}

// Static Constants (N) -------------------------------------------------------

#define NPOLY ((size_t) 1)
#define NBEAD ((size_t) 32)
#define NDIM ((size_t) 1)
#define NSUBDIV ((size_t) 16)

#define NTHRM ((size_t) 1 << 14)
#define NPROD ((size_t) 1 << 18)
#define NTHRMREC ((size_t) 1 << 4)
#define NPRODREC ((size_t) 1 << 8)

static void not_main(void) {
  err_reset();

  struct napkin* const napkin = napkin_alloc(NPOLY, NBEAD, NDIM, NSUBDIV,
      NTHRM, NPROD, NTHRMREC, NPRODREC);
  if (napkin == NULL)
    err_abort(napkin_alloc);

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

  napkin->periodic = true;
  napkin->L = 10.0;
  double const T = 0.1;
  napkin->beta = 1.0 / T;
  napkin->tau = napkin->beta / napkin->nmemb.bead;

  double const q = omega / 2.0;
  double const E = napkin->nmemb.dim * q / tanh(q * napkin->beta);
  (void) printf("Expected for QHO: E = %f\n", E);

  napkin->Vint = lj612b;
  napkin->Vend = zerob;
  napkin->Vext = harmb;

  conf_circlatt(napkin);
  force_close(napkin);
  // force_open(napkin);
  mass_const(napkin, 1.0);

  save_length(napkin);

  save_subdivisions(napkin);
  save_potential(napkin);

  FILE* const driftfp = fopen("qho-drift.data", "w");
  if (driftfp == NULL)
    err_abort(fopen);

  FILE* const tdefp = fopen("qho-tde.data", "w");
  if (tdefp == NULL)
    err_abort(fopen);

  // Just to prevent empty data files...
  disp_drift(napkin, driftfp);
  (void) fflush(driftfp);
  (void) fprintf(tdefp, "%zu %f %f %f\n", (size_t) 0, E, E, 0.0);
  (void) fflush(tdefp);

  static double (* const workers[])(struct napkin*) = {work_ss, work_comd};

  for (size_t istep = 0; istep < napkin->nstep.thrm + napkin->nstep.prod; ++istep) {
    choose(napkin,
        workers[ran_index(napkin->rng, sizeof workers / sizeof *workers)]
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

    if (istep < napkin->nstep.thrm) {
      if (napkin->nstep.prodrec * napkin->istep.thrm >
          napkin->nstep.prod * napkin->istep.thrmrec) {
        disp_drift(napkin, driftfp);

        ++napkin->istep.thrmrec;
      }

      ++napkin->istep.thrm;
    } else {
      est_accum(napkin->tde, est_tde(napkin));

      if (napkin->nstep.prodrec * napkin->istep.prod > napkin->nstep.prod * napkin->istep.prodrec) {
        disp_drift(napkin, driftfp);
        disp_tde(napkin, tdefp);

        ++napkin->istep.prodrec;
      }

      ++napkin->istep.prod;
    }
  }

  if (fclose(tdefp) == EOF)
    err_abort(fclose);

  if (fclose(driftfp) == EOF)
    err_abort(fclose);

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
// TODO Consider producing a histogram (half done).
// TODO Fix V/K scaling and offsets (done).
// TODO Abstract generic calculations for qho and He-4 (half done).
// TODO Figure out the correlation length for `sem` (half done).
// TODO Find home for lost souls (half done).
// TODO Consider different masses for different polymers (done).
// TODO PIGS time!

int main(void) {
  not_main();

  return EXIT_SUCCESS;
}
