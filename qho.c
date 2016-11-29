#include "err.h"
#include "fp.h"
#include "phys.h"
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

// Static Constants (N) -------------------------------------------------------

#define NPOLY ((size_t) 1)
#define NBEAD ((size_t) 64)
#define NDIM ((size_t) 1)
#define NSUBDIV ((size_t) 16)

#define NTSTEP ((size_t) 1 << 14)
#define NPSTEP ((size_t) 1 << 18)
#define NTRSTEP ((size_t) 1 << 4)
#define NPRSTEP ((size_t) 1 << 8)

// These assertions guarantee that adding or multiplying two steps or indices
// never wraps around, which makes their manipulation easier.
static_assert(NDIM * NPOLY * NBEAD <= SQRT_SIZE_MAX, "too many things");
static_assert(NTSTEP <= SQRT_SIZE_MAX, "too many thermalization steps");
static_assert(NPSTEP <= SQRT_SIZE_MAX, "too many production steps");
static_assert(NTRSTEP <= NTSTEP, "too many thermalization recording steps");
static_assert(NPRSTEP <= NPSTEP, "too many production recording steps");

// Types () -------------------------------------------------------------------

struct par {
  size_t accepted;
  size_t rejected;
  double dx;
};

struct estacc {
  size_t N; // Sample accumulator.
  double M1; // First moment accumulator.
  double M2; // Second moment accumulator.
};

struct estsum {
  double mean; // Mean.
  double var; // Variance.
  double sd; // Standard deviation.
  double sem; // Standard error of the mean.
};

struct step {
  size_t therm; // Thermalization steps.
  size_t prod; // Production steps.
  size_t thermrec; // Recording steps during thermalization.
  size_t prodrec; // Recording steps during production.
};

struct bead {
  double* d;
};

struct poly {
  double m;
  size_t from;
  size_t to;
  struct bead* r;
};

struct napkin {
  gsl_rng* rng; // State of the random number generator.
  struct step n; // Eventually...
  struct step i;
  double (* Vint)(struct napkin const*, struct bead, struct bead); // Potential between beads.
  double (* Vend)(struct napkin const*, struct bead, struct bead); // Extra potential between ends.
  double (* Vext)(struct napkin const*, struct bead); // External potential.
  double L; // Lattice constant.
  double beta; // Inverse temperature sometimes.
  double lambda; // De Broglie thing.
  double tau; // Imaginary time step maybe.
  struct {
    struct {
      size_t accepted;
      size_t rejected;
      double dx;
    } ssm; // Single slice move.
    struct {
      size_t accepted;
      size_t rejected;
      double dx;
    } comd; // Center of mass displacement.
  } params; // Optimization parameters of each move.
  struct {
    struct {
      size_t ipoly;
      size_t ibead;
      struct bead r;
    } ssm;
    struct {
      size_t ipoly;
      struct poly R;
    } comd;
  } history; // Undo history of the latest move.
  void (* accept)(struct napkin*); // Procedures for processing the latest move.
  void (* reject)(struct napkin*);
  void (* adjust)(struct napkin*);
  struct {
    struct estacc tde; // Thermodynamic energy estimator.
    struct estacc clp[NDIM]; // Central link position.
  } acc; // Estimator accumulators.
  size_t npoly;
  size_t nbead;
  size_t ndim;
  struct poly* R;
};

// Lost Souls () --------------------------------------------------------------

static size_t ran_index(gsl_rng* const rng, size_t const n) {
  return (size_t) gsl_rng_uniform_int(rng, n);
}

// Vector Math (bead) ---------------------------------------------------------

// Minimal norm squared from origin.
static double bead_norm2(struct napkin const* const napkin,
    struct bead const r) {
  double s = 0.0;

  for (size_t idim = 0; idim < NDIM; ++idim)
    s += gsl_pow_2(fp_wrap(r.d[idim], napkin->L));

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

  for (size_t idim = 0; idim < NDIM; ++idim)
    s += gsl_pow_2(fp_wrap(r1.d[idim] - r0.d[idim], napkin->L));

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
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    napkin->R[ipoly].m = m;
}

// Configuration Forcing (force) ----------------------------------------------

// Close the current configuration.
static void force_close(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    napkin->R[ipoly].from = ipoly;
    napkin->R[ipoly].to = ipoly;
  }
}

// Open the current configuration.
static void force_open(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    napkin->R[ipoly].from = SIZE_MAX;
    napkin->R[ipoly].to = SIZE_MAX;
  }
}

// Cycle the current configuration.
static void force_cycle(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    napkin->R[ipoly].from = size_uwrap(ipoly + NPOLY - 1, NPOLY);
    napkin->R[ipoly].to = size_uwrap(ipoly + 1, NPOLY);
  }
}

// Initial Configurations (conf) ----------------------------------------------

// Random initial configuration.
static void conf_rand(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      for (size_t idim = 0; idim < NDIM; ++idim)
        napkin->R[ipoly].r[ibead].d[idim] =
          gsl_ran_flat(napkin->rng, 0.0, napkin->L);
}

// Random point initial configuration.
static void conf_randpt(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t idim = 0; idim < NDIM; ++idim) {
      double const x = gsl_ran_flat(napkin->rng, 0, napkin->L);

      for (size_t ibead = 0; ibead < NBEAD; ++ibead)
        napkin->R[ipoly].r[ibead].d[idim] = x;
    }
}

// Random lattice initial configuration.
static void conf_randlatt(struct napkin* const napkin) {
  double const w = ceil(pow(NPOLY, 1.0 / NDIM));
  double const v = napkin->L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
      size_t i = ipoly;

      for (size_t idim = 0; idim < NDIM; ++idim) {
        size_div_t const z = size_div(i, n);

        napkin->R[ipoly].r[ibead].d[idim] =
          (double) z.rem * v + gsl_ran_flat(napkin->rng, 0, v);
        i = z.quot;
      }
    }
}

// Point lattice initial configuration.
static void conf_ptlatt(struct napkin* const napkin) {
  double const w = ceil(pow(NPOLY, 1.0 / NDIM));
  double const v = napkin->L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
      size_t i = ipoly;

      for (size_t idim = 0; idim < NDIM; ++idim) {
        size_div_t const z = size_div(i, n);

        napkin->R[ipoly].r[ibead].d[idim] = (double) z.rem * v;
        i = z.quot;
      }
    }
}

// Circle lattice initial configuration.
// The equations are based on the theory of Lissajous knots.
static void conf_circlatt(struct napkin* const napkin) {
  double const w = ceil(pow(NPOLY, 1.0 / NDIM));
  double const v = napkin->L / w;
  size_t const n = (size_t) w;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
      double const t = (double) ibead / NBEAD;

      size_t i = ipoly;

      for (size_t idim = 0; idim < NDIM; ++idim) {
        double const phi = (double) idim / NDIM;
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
  size_t const ipoly = napkin->history.ssm.ipoly;
  size_t const ibead = napkin->history.ssm.ibead;

  for (size_t idim = 0; idim < NDIM; ++idim)
    napkin->R[ipoly].r[ibead].d[idim] = napkin->history.ssm.r.d[idim];

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

  napkin->history.ssm.ipoly = ipoly;
  napkin->history.ssm.ibead = ibead;

  for (size_t idim = 0; idim < NDIM; ++idim) {
    napkin->history.ssm.r.d[idim] = napkin->R[ipoly].r[ibead].d[idim];

    napkin->R[ipoly].r[ibead].d[idim] = fp_uwrap(
        napkin->R[ipoly].r[ibead].d[idim] +
        napkin->params.ssm.dx * gsl_ran_flat(napkin->rng, -1.0, 1.0),
        napkin->L);
  }
}

static void move_accept_comd(struct napkin* const napkin) {
  ++napkin->params.comd.accepted;
}

static void move_reject_comd(struct napkin* const napkin) {
  size_t const ipoly = napkin->history.comd.ipoly;

  for (size_t ibead = 0; ibead < NBEAD; ++ibead)
    for (size_t idim = 0; idim < NDIM; ++idim)
      napkin->R[ipoly].r[ibead].d[idim] =
        napkin->history.comd.R.r[ibead].d[idim];

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

  napkin->history.comd.ipoly = ipoly;

  for (size_t idim = 0; idim < NDIM; ++idim) {
    double const x = napkin->params.comd.dx *
      gsl_ran_flat(napkin->rng, -1.0, 1.0);

    for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
      napkin->history.comd.R.r[ibead].d[idim] =
        napkin->R[ipoly].r[ibead].d[idim];

      napkin->R[ipoly].r[ibead].d[idim] = fp_uwrap(
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
  if (ibead == 0) {
    size_t const jpoly = napkin->R[ipoly].from;
    size_t const jbead = NBEAD - 1;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return bead_dist2(napkin,
          napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead - 1;

    return bead_dist2(napkin,
        napkin->R[ipoly].r[ibead], napkin->R[ipoly].r[jbead]);
  }
}

static double K_polybead_fw(struct napkin const* const napkin,
    size_t const ipoly, size_t const ibead) {
  if (ibead == NBEAD - 1) {
    size_t const jpoly = napkin->R[ipoly].to;
    size_t const jbead = 0;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return bead_dist2(napkin,
          napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead + 1;

    return bead_dist2(napkin,
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

  for (size_t ibead = 0; ibead < NBEAD; ++ibead)
    K += K_polybead_fw(napkin, ipoly, ibead);

  return K;
}

static double K_total(struct napkin const* const napkin) {
  double K = 0.0;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    K += K_poly(napkin, ipoly);

  return K;
}

static double Vint_bead(struct napkin const* const napkin,
    size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < NPOLY; ++jpoly)
      V += napkin->Vint(napkin,
          napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[ibead]);

  return V;
}

static double Vend_bead(struct napkin const* const napkin,
    size_t const ibead) {
  double V = 0.0;

  if (ibead == 0) {
    for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
      if (napkin->R[ipoly].from != SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < NPOLY; ++jpoly)
          if (napkin->R[jpoly].from != SIZE_MAX)
            V += napkin->Vend(napkin,
                  napkin->R[ipoly].r[ibead], napkin->R[jpoly].r[ibead]);
  } else if (ibead == NBEAD - 1)
    for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
      if (napkin->R[ipoly].to != SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < NPOLY; ++jpoly)
          if (napkin->R[jpoly].to != SIZE_MAX)
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

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly)
    V += Vext_polybead(napkin, ipoly, ibead);

  return V;
}

static double V_total(struct napkin const* const napkin) {
  double V = 0.0;

  for (size_t ibead = 0; ibead < NBEAD; ++ibead)
    V += V_bead(napkin, ibead);

  return V;
}

// Action Stuff (S) -----------------------------------------------------------

static double SV(struct napkin const* const napkin, double const V) {
  return napkin->tau * V;
}

static double SK(struct napkin const* const napkin, double const K) {
  return K / (4.0 * napkin->lambda * napkin->tau);
}

// Estimator Statistics (est) -------------------------------------------------

static struct estacc est_empty(void) {
  struct estacc acc;

  acc.N = 0;
  acc.M1 = 0.0;
  acc.M2 = 0.0;

  return acc;
}

static struct estacc est_accum(struct estacc acc, double const x) {
  ++acc.N;
  double const dx = x - acc.M1;
  acc.M1 += dx / (double) acc.N;
  acc.M2 += dx * (x - acc.M1);

  return acc;
}

static struct estsum est_summ(struct estacc const acc) {
  struct estsum sum;

  switch (acc.N) {
    case 0:
      sum.mean = NAN;
      sum.var = NAN;
      sum.sd = NAN;
      sum.sem = NAN;
      break;
    case 1:
      sum.mean = acc.M1;
      sum.var = 0.0;
      sum.sd = 0.0;
      sum.sem = 0.0;
      break;
    default:
      sum.mean = acc.M1;
      sum.var = acc.M2 / (double) (acc.N - 1);
      sum.sd = sqrt(sum.var);
      sum.sem = sum.sd / sqrt((double) acc.N);
  }

  return sum;
}

// Estimators (est) -----------------------------------------------------------

// \langle E_T\rangle & = \frac 1{M \tau} \sum_{k = 1}^M
// \Bigl\langle\frac{d N} 2 - \frac{|R_{k + 1 \bmod M} - R_k|^2}{4 \lambda \tau} + \tau V(R_k)\Bigr\rangle
static double est_tde(struct napkin const* const napkin) {
  double const KC = NDIM * NPOLY * NBEAD / 2.0;
  double const K = SK(napkin, K_total(napkin));
  double const V = SV(napkin, V_total(napkin));

  return (KC - K + V) / (NPOLY * NBEAD * napkin->tau);
}

// Workers (work) -------------------------------------------------------------

static double work_ss(struct napkin* const napkin) {
  size_t const ipoly = ran_index(napkin->rng, NPOLY);
  size_t const ibead = ran_index(napkin->rng, NBEAD);

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
  size_t const ipoly = ran_index(napkin->rng, NPOLY);

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

  for (size_t idim = 0; idim < NDIM; ++idim)
    (void) fprintf(fp, " %f", napkin->R[ipoly].r[ibead].d[idim]);

  (void) fprintf(fp, "\n");
}

static void disp_poly(struct napkin const* const napkin,
    FILE* const fp) {
  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      disp_bead(napkin, fp, ibead, ipoly, ibead);

    if (napkin->R[ipoly].to != SIZE_MAX)
      disp_bead(napkin, fp, NBEAD, napkin->R[ipoly].to, 0);

    (void) fprintf(fp, "\n");
  }
}

static void disp_tde(struct napkin const* const napkin,
    FILE* const fp) {
  struct estsum const sum = est_summ(napkin->acc.tde);

  (void) fprintf(fp, "%zu %f %f %f\n",
      // Correlations!
      napkin->i.prod, est_tde(napkin), sum.mean, sum.sem * NBEAD);
}

static void disp_drift(struct napkin const* const napkin,
    FILE* const fp) {
  (void) fprintf(fp, "%zu %f %f\n",
      napkin->i.therm + napkin->i.prod,
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
  if (snprintf(str, BUFSIZ, "qho-ensemble-%zud.data", NDIM) < 0)
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

  (void) fprintf(fp, "%zu\n", NSUBDIV);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

static void save_potential(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-potential-%zud.data", NDIM) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  double const v = napkin->L / NSUBDIV;

  double* const d = malloc(NDIM * sizeof *d);
  if (d == NULL)
    err_abort(malloc);

  struct bead r = {.d = d};

  for (size_t ipt = 0; ipt < size_pow(NSUBDIV + 1, NDIM); ++ipt) {
    size_t i = ipt;

    for (size_t idim = 0; idim < NDIM; ++idim) {
      size_div_t const z = size_div(i, NSUBDIV + 1);

      r.d[idim] = (double) z.rem * v;
      i = z.quot;
    }

    for (size_t idim = 0; idim < NDIM; ++idim)
      (void) fprintf(fp, "%f ", r.d[idim]);

    (void) fprintf(fp, "%f\n", napkin->Vext(napkin, r));
  }

  free(d);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

// Other Crap () --------------------------------------------------------------

static void status_line(struct napkin const* const napkin) {
  struct estsum const sum = est_summ(napkin->acc.tde);

  (void) printf("i / N = (%zu + %zu) / (%zu + %zu) = %zu %%\n"
      "E = %f +- %f\n",
      napkin->i.therm, napkin->i.prod, napkin->n.therm, napkin->n.prod,
      100 *
      (napkin->i.therm + napkin->i.prod) / (napkin->n.therm + napkin->n.prod),
      // Correlations!
      sum.mean, sum.sem * NBEAD);
}

static void status_file(struct napkin const* const napkin) {
  save_esl(napkin);
}

// More Crap () ---------------------------------------------------------------

static double lj612b(struct napkin const* const napkin,
    struct bead const r0, struct bead const r1) {
  double const epsilon = 1.0;
  double const sigma = 1.0;

  double const sigmar2 = gsl_pow_2(sigma) / bead_dist2(napkin, r0, r1);

  return 4.0 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

static double harmb(struct napkin const* const napkin,
    struct bead const r) {
  double const m = 1.0;
  double const omega = 1.0;

  return m * gsl_pow_2(omega) * bead_norm2(napkin, r) / 2.0;
}

static double zerob(__attribute__ ((__unused__))
    struct napkin const* const napkin,
    __attribute__ ((__unused__)) struct bead const r0,
    __attribute__ ((__unused__)) struct bead const r1) {
  return 0.0;
}

static void not_main(void) {
  err_reset();

  struct napkin* const napkin = malloc(sizeof *napkin);
  if (napkin == NULL)
    err_abort(malloc);

  gsl_rng_type const* const t = gsl_rng_env_setup();
  napkin->rng = gsl_rng_alloc(t);
  if (napkin->rng == NULL)
    err_abort(gsl_rng_alloc);

  napkin->R = malloc(NPOLY * sizeof *napkin->R);
  if (napkin->R == NULL)
    err_abort(malloc);

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    napkin->R[ipoly].r = malloc(NBEAD * sizeof *napkin->R[ipoly].r);
    if (napkin->R[ipoly].r == NULL)
      err_abort(malloc);

    for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
      napkin->R[ipoly].r[ibead].d =
        malloc(NDIM * sizeof *napkin->R[ipoly].r[ibead].d);
      if (napkin->R[ipoly].r[ibead].d == NULL)
        err_abort(malloc);
    }
  }

  napkin->history.ssm.r.d = malloc(NDIM * sizeof *napkin->history.ssm.r.d);
  if (napkin->history.ssm.r.d == NULL)
    err_abort(malloc);

  napkin->history.comd.R.r = malloc(NBEAD * sizeof *napkin->history.comd.R.r);
  if (napkin->history.comd.R.r == NULL)
    err_abort(malloc);

  for (size_t ibead = 0; ibead < NBEAD; ++ibead) {
    napkin->history.comd.R.r[ibead].d =
      malloc(NDIM * sizeof *napkin->history.comd.R.r[ibead].d);
    if (napkin->history.comd.R.r[ibead].d == NULL)
      err_abort(malloc);
  }

  int const sigs[] = {SIGUSR1, SIGUSR2};
  if (sigs_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX)
    err_abort(sigs_register);

  napkin->params.ssm.accepted = 0;
  napkin->params.ssm.rejected = 0;
  napkin->params.ssm.dx = 1;

  napkin->params.comd.accepted = 0;
  napkin->params.comd.rejected = 0;
  napkin->params.comd.dx = 1;

  double const T = 0.1;
  double const m = 1.0;
  napkin->beta = 1.0 / (kB * T);
  napkin->tau = napkin->beta / NBEAD;
  napkin->L = 10.0;

  napkin->lambda = gsl_pow_2(hbar) / (2.0 * m); // TODO Per m.

  double const omega = 1.0;
  double const q = hbar * omega / 2.0;
  double const E = q / tanh(q * napkin->beta);
  (void) printf("expect in 1d: E = %f\n", E);

  napkin->Vint = lj612b;
  napkin->Vend = zerob;
  napkin->Vext = harmb;

  napkin->acc.tde = est_empty();
  for (size_t idim = 0; idim < NDIM; ++idim)
    napkin->acc.clp[idim] = est_empty();

  conf_circlatt(napkin);
  force_close(napkin);
  // force_open(napkin);
  mass_const(napkin, m);

  save_length(napkin);

  save_subdivisions(napkin);
  save_potential(napkin);

  FILE* const driftfp = fopen("qho-drift.data", "w");
  if (driftfp == NULL)
    err_abort(fopen);

  FILE* const tdefp = fopen("qho-tde.data", "w");
  if (tdefp == NULL)
    err_abort(fopen);

  napkin->n.therm = NTSTEP;
  napkin->n.prod = NPSTEP;
  napkin->n.thermrec = NTRSTEP;
  napkin->n.prodrec = NPRSTEP;

  napkin->i.therm = 0;
  napkin->i.prod = 0;
  napkin->i.thermrec = 0;
  napkin->i.prodrec = 0;

  // Just to prevent empty data files...
  disp_drift(napkin, driftfp);
  (void) fflush(driftfp);
  (void) fprintf(tdefp, "%zu %f %f %f\n", (size_t) 0, E, E, 0.0);
  (void) fflush(tdefp);

  static double (* const workers[])(struct napkin*) = {work_ss, work_comd};

  for (size_t istep = 0; istep < napkin->n.therm + napkin->n.prod; ++istep) {
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

    if (istep < napkin->n.therm) {
      if (napkin->n.prodrec * napkin->i.therm > napkin->n.prod * napkin->i.thermrec) {
        disp_drift(napkin, driftfp);

        ++napkin->i.thermrec;
      }

      ++napkin->i.therm;
    } else {
      napkin->acc.tde = est_accum(napkin->acc.tde, est_tde(napkin));

      for (size_t idim = 0; idim < NDIM; ++idim)
        napkin->acc.clp[idim] = est_accum(napkin->acc.clp[idim],
            fp_wrap(napkin->R[0].r[NBEAD / 2].d[idim], napkin->L));

      if (napkin->n.prodrec * napkin->i.prod > napkin->n.prod * napkin->i.prodrec) {
        disp_drift(napkin, driftfp);
        disp_tde(napkin, tdefp);

        ++napkin->i.prodrec;
      }

      ++napkin->i.prod;
    }
  }

  // Bad.
  for (size_t idim = 0; idim < NDIM; ++idim) {
    struct estsum const sum = est_summ(napkin->acc.clp[idim]);

    (void) printf("should be zero (%zu / %zu): %f +- %f\n",
        idim + 1, NDIM, sum.mean, sum.sem);
  }

  if (fclose(tdefp) == EOF)
    err_abort(fclose);

  if (fclose(driftfp) == EOF)
    err_abort(fclose);

  status_line(napkin);
  status_file(napkin);

  for (size_t ibead = 0; ibead < NBEAD; ++ibead)
    free(napkin->history.comd.R.r[ibead].d);

  free(napkin->history.comd.R.r);

  free(napkin->history.ssm.r.d);

  for (size_t ipoly = 0; ipoly < NPOLY; ++ipoly) {
    for (size_t ibead = 0; ibead < NBEAD; ++ibead)
      free(napkin->R[ipoly].r[ibead].d);

    free(napkin->R[ipoly].r);
  }

  free(napkin->R);

  gsl_rng_free(napkin->rng);

  free(napkin);
}

// Main and Procrastination () ------------------------------------------------

// TODO Try other Marsaglia's rngs (UNI or VNI had no observable effect, done).
// TODO Split data files by dimension (use `snprintf`, done).
// TODO Abstract mean/var (done).
// TODO Meditate about the anisotropy and symmetry of potentials (done).
// TODO Deal with the disparity `/ 2` in K_total and V_total (done).
// TODO Make periodicity conditional (half done).
// TODO Consider producing a histogram (half done).
// TODO Fix V/K scaling and offsets (done).
// TODO Abstract generic calculations for qho and He-4 (half done).
// TODO Figure out the correlation length for `sem` (half done).
// TODO Find home for lost souls (half done).
// TODO Consider different masses for different polymers (half done).
// TODO PIGS time!

int main(void) {
  not_main();

  return EXIT_SUCCESS;
}
