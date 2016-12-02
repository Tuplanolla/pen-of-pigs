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
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Types () -------------------------------------------------------------------

// This structure contains the number of thermalization and production steps
// in addition to the number of recording steps for both.
struct step {
  size_t thrm;
  size_t prod;
  size_t thrmrec;
  size_t prodrec;
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

struct ensem {
  struct poly* R;
  double (* Vint)(struct ensem const*, struct bead, struct bead);
  double (* Vend)(struct ensem const*, struct bead, struct bead);
  double (* Vext)(struct ensem const*, struct bead);
  struct memb nmemb;
  struct step nstep;
  struct step istep;
  double L;
  double tau;
  bool periodic;
};

struct napkin {
  gsl_rng* rng;
  struct stats* tde;
  void (* accept)(struct napkin*);
  void (* reject)(struct napkin*);
  void (* adjust)(struct napkin*);
  struct ensem ensem;
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
    } comd;
  } params;
};

// Adjustment Surfaces () -----------------------------------------------------

static double surface(double const a, double const r) {
  return 0.5 + a / (a + r);
}

__attribute__ ((__const__))
static double another_surface(double const a, double const r) {
  return 1.0 - exp(-a) + exp(-r);
}

// Vector Math (bead) ---------------------------------------------------------

// Minimal norm squared from origin.
__attribute__ ((__nonnull__, __pure__))
static double bead_norm2(struct ensem const* const ensem,
    struct bead const r) {
  double s = 0.0;

  double (* const f)(double, double) =
    ensem->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ensem->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r.d[idim], ensem->L));

  return s;
}

// Minimal norm from origin.
__attribute__ ((__nonnull__, __pure__))
static double bead_norm(struct ensem const* const ensem,
    struct bead const r) {
  return sqrt(bead_norm2(ensem, r));
}

// Minimal distance squared between two points.
__attribute__ ((__nonnull__, __pure__))
static double bead_dist2(struct ensem const* const ensem,
    struct bead const r0, struct bead const r1) {
  double s = 0.0;

  double (* const f)(double, double) =
    ensem->periodic ? fp_wrap : fp_constant;

  for (size_t idim = 0; idim < ensem->nmemb.dim; ++idim)
    s += gsl_pow_2(f(r1.d[idim] - r0.d[idim], ensem->L));

  return s;
}

// Minimal distance between two points.
__attribute__ ((__nonnull__, __pure__))
static double bead_dist(struct ensem const* const ensem,
    struct bead const r0, struct bead const r1) {
  return sqrt(bead_dist2(ensem, r0, r1));
}

// Kinetic and Potential Energies (K, V) --------------------------------------

__attribute__ ((__nonnull__, __pure__))
static double K_polybead_bw(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  double const M = ensem->R[ipoly].m / (2.0 * gsl_pow_2(ensem->tau));

  if (ibead == 0) {
    size_t const jpoly = ensem->R[ipoly].from;
    size_t const jbead = ensem->nmemb.bead - 1;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return M * bead_dist2(ensem,
          ensem->R[ipoly].r[ibead], ensem->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead - 1;

    return M * bead_dist2(ensem,
        ensem->R[ipoly].r[ibead], ensem->R[ipoly].r[jbead]);
  }
}

// \lambda = \hbar / 2 m, K = x^2 / 4 \lambda \tau^2 = M x^2
// M = m / 2 \hbar \tau^2
__attribute__ ((__nonnull__, __pure__))
static double K_polybead_fw(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  double const M = ensem->R[ipoly].m / (2.0 * gsl_pow_2(ensem->tau));

  if (ibead == ensem->nmemb.bead - 1) {
    size_t const jpoly = ensem->R[ipoly].to;
    size_t const jbead = 0;

    if (jpoly == SIZE_MAX)
      return 0.0;
    else
      return M * bead_dist2(ensem,
          ensem->R[ipoly].r[ibead], ensem->R[jpoly].r[jbead]);
  } else {
    size_t const jbead = ibead + 1;

    return M * bead_dist2(ensem,
        ensem->R[ipoly].r[ibead], ensem->R[ipoly].r[jbead]);
  }
}

__attribute__ ((__nonnull__, __pure__))
static double K_polybead(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  return K_polybead_bw(ensem, ipoly, ibead) +
    K_polybead_fw(ensem, ipoly, ibead);
}

__attribute__ ((__nonnull__, __pure__))
static double K_poly(struct ensem const* const ensem,
    size_t const ipoly) {
  double K = 0.0;

  for (size_t ibead = 0; ibead < ensem->nmemb.bead; ++ibead)
    K += K_polybead_fw(ensem, ipoly, ibead);

  return K;
}

__attribute__ ((__nonnull__, __pure__))
static double K_total(struct ensem const* const ensem) {
  double K = 0.0;

  for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
    K += K_poly(ensem, ipoly);

  return K;
}

__attribute__ ((__nonnull__, __pure__))
static double Vint_bead(struct ensem const* const ensem,
    size_t const ibead) {
  double V = 0.0;

  for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
    for (size_t jpoly = ipoly + 1; jpoly < ensem->nmemb.poly; ++jpoly)
      V += ensem->Vint(ensem,
          ensem->R[ipoly].r[ibead], ensem->R[jpoly].r[ibead]);

  return V;
}

__attribute__ ((__nonnull__, __pure__))
static double Vend_bead(struct ensem const* const ensem,
    size_t const ibead) {
  double V = 0.0;

  if (ibead == 0) {
    for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
      if (ensem->R[ipoly].from == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < ensem->nmemb.poly; ++jpoly)
          if (ensem->R[jpoly].from == SIZE_MAX)
            V += ensem->Vend(ensem,
                ensem->R[ipoly].r[ibead], ensem->R[jpoly].r[ibead]);
  } else if (ibead == ensem->nmemb.bead - 1)
    for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
      if (ensem->R[ipoly].to == SIZE_MAX)
        for (size_t jpoly = ipoly + 1; jpoly < ensem->nmemb.poly; ++jpoly)
          if (ensem->R[jpoly].to == SIZE_MAX)
            V += ensem->Vend(ensem,
                ensem->R[ipoly].r[ibead], ensem->R[jpoly].r[ibead]);

  return V;
}

__attribute__ ((__nonnull__, __pure__))
static double Vext_polybead(struct ensem const* const ensem,
    size_t const ipoly, size_t const ibead) {
  return ensem->Vext(ensem, ensem->R[ipoly].r[ibead]);
}

__attribute__ ((__nonnull__, __pure__))
static double V_bead(struct ensem const* const ensem,
    size_t const ibead) {
  double V = Vint_bead(ensem, ibead) + Vend_bead(ensem, ibead);

  for (size_t ipoly = 0; ipoly < ensem->nmemb.poly; ++ipoly)
    V += Vext_polybead(ensem, ipoly, ibead);

  return V;
}

__attribute__ ((__nonnull__, __pure__))
static double V_total(struct ensem const* const ensem) {
  double V = 0.0;

  for (size_t ibead = 0; ibead < ensem->nmemb.bead; ++ibead)
    V += V_bead(ensem, ibead);

  return V;
}

// Estimators (est) -----------------------------------------------------------

// Thermodynamic energy estimator.
// \langle E_T\rangle & = \frac 1{M \tau} \sum_{k = 1}^M
// \Bigl\langle\frac{d N} 2 - \frac{|R_{k + 1 \bmod M} - R_k|^2}{4 \lambda \tau} + \tau V(R_k)\Bigr\rangle
__attribute__ ((__nonnull__, __pure__))
static double est_tde(struct ensem const* const ensem) {
  double const E =
    (double) (ensem->nmemb.dim * ensem->nmemb.poly) / (2.0 * ensem->tau);
  double const K =
    K_total(ensem) / (double) (ensem->nmemb.poly * ensem->nmemb.bead);
  double const V =
    V_total(ensem) / (double) (ensem->nmemb.poly * ensem->nmemb.bead);

  return E - K + V;
}

// Mass Configurations (Rm) ---------------------------------------------------

// Equal mass for every particle.
__attribute__ ((__nonnull__))
static void Rm_const(struct napkin* const napkin, double const m) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    napkin->ensem.R[ipoly].m = m;
}

// End Configurations (Rend) --------------------------------------------------

// Close the current configuration.
__attribute__ ((__nonnull__))
static void Rend_close(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].from = ipoly;
    napkin->ensem.R[ipoly].to = ipoly;
  }
}

// Open the current configuration.
__attribute__ ((__nonnull__))
static void Rend_open(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].from = SIZE_MAX;
    napkin->ensem.R[ipoly].to = SIZE_MAX;
  }
}

// Cycle the current configuration.
__attribute__ ((__nonnull__))
static void Rend_cycle(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    napkin->ensem.R[ipoly].from = size_uwrap_dec(ipoly, napkin->ensem.nmemb.poly);
    napkin->ensem.R[ipoly].to = size_uwrap_inc(ipoly, napkin->ensem.nmemb.poly);
  }
}

// Bead Configurations (Rr) ---------------------------------------------------

// Random initial configuration.
__attribute__ ((__nonnull__))
static void Rr_rand(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
        napkin->ensem.R[ipoly].r[ibead].d[idim] =
          ran_uopen(napkin->rng, napkin->ensem.L);
}

// Random point initial configuration.
__attribute__ ((__nonnull__))
static void Rr_randpt(struct napkin* const napkin) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
      double const x = ran_uopen(napkin->rng, napkin->ensem.L);

      for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
        napkin->ensem.R[ipoly].r[ibead].d[idim] = x;
    }
}

// Random lattice initial configuration.
__attribute__ ((__nonnull__))
static void Rr_randlatt(struct napkin* const napkin) {
  size_t const n = size_cirt(napkin->ensem.nmemb.poly, napkin->ensem.nmemb.dim);
  double const v = napkin->ensem.L / (double) n;

  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
      size_t i = ipoly;

      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
        size_div_t const z = size_div(i, n);

        napkin->ensem.R[ipoly].r[ibead].d[idim] =
          (double) z.rem * v + ran_uopen(napkin->rng, v);
        i = z.quot;
      }
    }
}

// Point lattice initial configuration.
__attribute__ ((__nonnull__))
static void Rr_ptlatt(struct napkin* const napkin) {
  size_t const n = size_cirt(napkin->ensem.nmemb.poly, napkin->ensem.nmemb.dim);
  double const v = napkin->ensem.L / (double) n;

  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
      size_t i = ipoly;

      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
        size_div_t const z = size_div(i, n);

        napkin->ensem.R[ipoly].r[ibead].d[idim] = (double) z.rem * v;
        i = z.quot;
      }
    }
}

// Circle lattice initial configuration.
// The equations are based on the theory of Lissajous knots.
__attribute__ ((__nonnull__))
static void Rr_circlatt(struct napkin* const napkin) {
  size_t const n = size_cirt(napkin->ensem.nmemb.poly, napkin->ensem.nmemb.dim);
  double const v = napkin->ensem.L / (double) n;

  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly)
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
      double const t = (double) ibead / (double) napkin->ensem.nmemb.bead;

      size_t i = ipoly;

      for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
        double const phi = (double) idim / (double) napkin->ensem.nmemb.dim;
        double const y = sin(M_2PI * (t + phi / 2.0)) / 2.0;

        size_div_t const z = size_div(i, n);
        napkin->ensem.R[ipoly].r[ibead].d[idim] =
          ((double) z.rem + (y + 1.0) / 2.0) * v;
        i = z.quot;
      }
    }
}

// Physical Moves (move) ------------------------------------------------------

// r <=> a: z >=< 1

__attribute__ ((__const__))
__attribute__ ((__nonnull__))
static void move_accept_ss(struct napkin* const napkin) {
  ++napkin->params.ssm.accepted;
}

__attribute__ ((__nonnull__))
static void move_reject_ss(struct napkin* const napkin) {
  size_t const ipoly = napkin->hist.ssm.ipoly;
  size_t const ibead = napkin->hist.ssm.ibead;

  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
    napkin->ensem.R[ipoly].r[ibead].d[idim] = napkin->hist.ssm.r.d[idim];

  ++napkin->params.ssm.rejected;
}

__attribute__ ((__nonnull__))
static void move_adjust_ss(struct napkin* const napkin) {
  napkin->params.ssm.dx = napkin->params.ssm.dx *
    surface((double) napkin->params.ssm.accepted,
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

  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
    napkin->hist.ssm.r.d[idim] = napkin->ensem.R[ipoly].r[ibead].d[idim];

    napkin->ensem.R[ipoly].r[ibead].d[idim] =
      f(napkin->ensem.R[ipoly].r[ibead].d[idim] +
        napkin->params.ssm.dx * ran_open(napkin->rng, 1.0), napkin->ensem.L);
  }
}

__attribute__ ((__nonnull__))
static void move_accept_comd(struct napkin* const napkin) {
  ++napkin->params.comd.accepted;
}

__attribute__ ((__nonnull__))
static void move_reject_comd(struct napkin* const napkin) {
  size_t const ipoly = napkin->hist.comd.ipoly;

  for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
      napkin->ensem.R[ipoly].r[ibead].d[idim] =
        napkin->hist.comd.R.r[ibead].d[idim];

  ++napkin->params.comd.rejected;
}

__attribute__ ((__nonnull__))
static void move_adjust_comd(struct napkin* const napkin) {
  napkin->params.comd.dx = napkin->params.comd.dx *
    another_surface((double) napkin->params.comd.accepted,
        (double) napkin->params.comd.rejected);
}

__attribute__ ((__nonnull__))
static void move_comd(struct napkin* const napkin,
    size_t const ipoly) {
  napkin->accept = move_accept_comd;
  napkin->reject = move_reject_comd;
  napkin->adjust = move_adjust_comd;

  napkin->hist.comd.ipoly = ipoly;

  double (* const f)(double, double) =
    napkin->ensem.periodic ? fp_uwrap : fp_constant;

  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
    double const x = napkin->params.comd.dx *
      ran_open(napkin->rng, 1.0);

    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
      napkin->hist.comd.R.r[ibead].d[idim] =
        napkin->ensem.R[ipoly].r[ibead].d[idim];

      napkin->ensem.R[ipoly].r[ibead].d[idim] =
        f(napkin->ensem.R[ipoly].r[ibead].d[idim] + x, napkin->ensem.L);
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

// Workers (work) -------------------------------------------------------------

__attribute__ ((__nonnull__))
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

__attribute__ ((__nonnull__))
static double work_comd(struct napkin* const napkin) {
  size_t const ipoly = ran_index(napkin->rng, napkin->ensem.nmemb.poly);

  double const V0 = V_total(&napkin->ensem);

  move_comd(napkin, ipoly);

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

// Printing into Data Files (disp) --------------------------------------------

__attribute__ ((__nonnull__))
static void disp_bead(struct napkin const* const napkin,
    FILE* const fp,
    size_t const iindex, size_t const ipoly, size_t const ibead) {
  (void) fprintf(fp, "%zu", iindex);

  for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
    (void) fprintf(fp, " %f", napkin->ensem.R[ipoly].r[ibead].d[idim]);

  (void) fprintf(fp, "\n");
}

__attribute__ ((__nonnull__))
static void disp_poly(struct napkin const* const napkin,
    FILE* const fp) {
  for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
    for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
      disp_bead(napkin, fp, ibead, ipoly, ibead);

    if (napkin->ensem.R[ipoly].to != SIZE_MAX)
      disp_bead(napkin, fp, napkin->ensem.nmemb.bead, napkin->ensem.R[ipoly].to, 0);

    (void) fprintf(fp, "\n");
  }
}

__attribute__ ((__nonnull__))
static void disp_tde(struct napkin const* const napkin,
    FILE* const fp) {
  (void) fprintf(fp, "%zu %f %f %f\n",
      napkin->ensem.istep.prod, est_tde(&napkin->ensem), stats_mean(napkin->tde),
      // Correlations!
      stats_sem(napkin->tde) * (double) napkin->ensem.nmemb.bead);
}

__attribute__ ((__nonnull__))
static void disp_drift(struct napkin const* const napkin,
    FILE* const fp) {
  (void) fprintf(fp, "%zu %f %f\n",
      napkin->ensem.istep.thrm + napkin->ensem.istep.prod,
      napkin->params.ssm.dx, napkin->params.comd.dx);
}

// Saving into Data Files (save)

__attribute__ ((__nonnull__))
static void save_length(struct napkin const* const napkin) {
  FILE* const fp = fopen("qho-length.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  (void) fprintf(fp, "%f\n", napkin->ensem.L);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_esl(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-ensemble-%zud.data", napkin->ensem.nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  disp_poly(napkin, fp);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_subdivisions(struct napkin const* const napkin) {
  FILE* const fp = fopen("qho-subdivisions.data", "w");
  if (fp == NULL)
    err_abort(fopen);

  (void) fprintf(fp, "%zu\n", napkin->ensem.nmemb.subdiv);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

__attribute__ ((__nonnull__))
static void save_potential(struct napkin const* const napkin) {
  char str[BUFSIZ]; // Always longer than 256 and safe!
  if (snprintf(str, BUFSIZ, "qho-potential-%zud.data", napkin->ensem.nmemb.dim) < 0)
    err_abort(snprintf);

  FILE* const fp = fopen(str, "w");
  if (fp == NULL)
    err_abort(fopen);

  double const v = napkin->ensem.L / (double) napkin->ensem.nmemb.subdiv;

  double* const d = malloc(napkin->ensem.nmemb.dim * sizeof *d);
  if (d == NULL)
    err_abort(malloc);

  struct bead r = {.d = d};

  for (size_t ipt = 0; ipt < size_pow(napkin->ensem.nmemb.subdiv + 1, napkin->ensem.nmemb.dim); ++ipt) {
    size_t i = ipt;

    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim) {
      size_div_t const z = size_div(i, napkin->ensem.nmemb.subdiv + 1);

      r.d[idim] = (double) z.rem * v;
      i = z.quot;
    }

    for (size_t idim = 0; idim < napkin->ensem.nmemb.dim; ++idim)
      (void) fprintf(fp, "%f ", r.d[idim]);

    (void) fprintf(fp, "%f\n", napkin->ensem.Vext(&napkin->ensem, r));
  }

  free(d);

  if (fclose(fp) == EOF)
    err_abort(fclose);
}

// Other Crap () --------------------------------------------------------------

__attribute__ ((__nonnull__))
static void status_line(struct napkin const* const napkin) {
  (void) printf("i / N = (%zu + %zu) / (%zu + %zu) = %zu %%\n"
      "E = %f +- %f\n",
      napkin->ensem.istep.thrm, napkin->ensem.istep.prod, napkin->ensem.nstep.thrm, napkin->ensem.nstep.prod,
      100 * (napkin->ensem.istep.thrm + napkin->ensem.istep.prod) /
      (napkin->ensem.nstep.thrm + napkin->ensem.nstep.prod),
      stats_mean(napkin->tde),
      // Correlations!
      stats_sem(napkin->tde) * (double) napkin->ensem.nmemb.bead);
}

__attribute__ ((__nonnull__))
static void status_file(struct napkin const* const napkin) {
  save_esl(napkin);
}

// Allocation (alloc) ---------------------------------------------------------

typedef void* (* alloc_proc)(size_t);

static void* alloc_err(__attribute__ ((__unused__)) size_t const n) {
  return NULL;
}

// Napkin Management (napkin) -------------------------------------------------

__attribute__ ((__nonnull__))
static void napkin_free(struct napkin* const napkin) {
  if (napkin != NULL) {
    if (napkin->hist.comd.R.r != NULL)
      for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
        free(napkin->hist.comd.R.r[ibead].d);

    free(napkin->hist.comd.R.r);

    free(napkin->hist.ssm.r.d);

    stats_free(napkin->tde);

    if (napkin->ensem.R != NULL)
      for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
        if (napkin->ensem.R[ipoly].r != NULL)
          for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead)
            free(napkin->ensem.R[ipoly].r[ibead].d);

        free(napkin->ensem.R[ipoly].r);
      }

    free(napkin->ensem.R);

    gsl_rng_free(napkin->rng);
  }

  free(napkin);
}

__attribute__ ((__malloc__))
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

  alloc_proc alloc = malloc;

  struct napkin* const napkin = alloc(sizeof *napkin);
  if (napkin == NULL)
    alloc = alloc_err;
  else {
    napkin->ensem.nmemb.poly = npoly;
    napkin->ensem.nmemb.bead = nbead;
    napkin->ensem.nmemb.dim = ndim;
    napkin->ensem.nmemb.subdiv = nsubdiv;

    napkin->ensem.nstep.thrm = nthrm;
    napkin->ensem.nstep.prod = nprod;
    napkin->ensem.nstep.thrmrec = nthrm;
    napkin->ensem.nstep.prodrec = nprodrec;

    napkin->ensem.istep.thrm = 0;
    napkin->ensem.istep.prod = 0;
    napkin->ensem.istep.thrmrec = 0;
    napkin->ensem.istep.prodrec = 0;

    napkin->rng = gsl_rng_alloc(gsl_rng_env_setup());
    if (napkin->rng == NULL)
      alloc = alloc_err;

    napkin->ensem.R = alloc(napkin->ensem.nmemb.poly * sizeof *napkin->ensem.R);
    if (napkin->ensem.R == NULL)
      alloc = alloc_err;
    else
      for (size_t ipoly = 0; ipoly < napkin->ensem.nmemb.poly; ++ipoly) {
        napkin->ensem.R[ipoly].r =
          alloc(napkin->ensem.nmemb.bead * sizeof *napkin->ensem.R[ipoly].r);
        if (napkin->ensem.R[ipoly].r == NULL)
          alloc = alloc_err;
        else
          for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
            napkin->ensem.R[ipoly].r[ibead].d =
              alloc(napkin->ensem.nmemb.dim * sizeof *napkin->ensem.R[ipoly].r[ibead].d);
            if (napkin->ensem.R[ipoly].r[ibead].d == NULL)
              alloc = alloc_err;
          }
      }

    napkin->hist.ssm.r.d =
      alloc(napkin->ensem.nmemb.dim * sizeof *napkin->hist.ssm.r.d);
    if (napkin->hist.ssm.r.d == NULL)
      alloc = alloc_err;

    napkin->hist.comd.R.r =
      alloc(napkin->ensem.nmemb.bead * sizeof *napkin->hist.comd.R.r);
    if (napkin->hist.comd.R.r == NULL)
      alloc = alloc_err;
    else
      for (size_t ibead = 0; ibead < napkin->ensem.nmemb.bead; ++ibead) {
        napkin->hist.comd.R.r[ibead].d =
          alloc(napkin->ensem.nmemb.dim * sizeof *napkin->hist.comd.R.r[ibead].d);
        if (napkin->hist.comd.R.r[ibead].d == NULL)
          alloc = alloc_err;
      }

    napkin->tde = stats_alloc();
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

// More Crap () ---------------------------------------------------------------

static double const epsilon = 1.0;
static double const sigma = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double V_lj612(struct ensem const* const ensem,
    struct bead const r0, struct bead const r1) {
  double const sigmar2 = gsl_pow_2(sigma) / bead_dist2(ensem, r0, r1);

  return 4.0 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

static double const omega = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double Vext_harm(struct ensem const* const ensem, struct bead const r) {
  return ensem->R[0].m * gsl_pow_2(omega) * bead_norm2(ensem, r) / 2.0;
}

__attribute__ ((__const__, __nonnull__))
static double V_zero(
    __attribute__ ((__unused__)) struct ensem const* const ensem,
    __attribute__ ((__unused__)) struct bead const r0,
    __attribute__ ((__unused__)) struct bead const r1) {
  return 0.0;
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

  napkin->ensem.periodic = true;
  napkin->ensem.L = 10.0;
  double const T = 0.1;
  double const beta = 1.0 / T;
  napkin->ensem.tau = beta / (double) napkin->ensem.nmemb.bead;

  double const q = omega / 2.0;
  double const E = (double) napkin->ensem.nmemb.dim * q / tanh(q * beta);
  (void) printf("Expected for QHO: E = %f\n", E);

  napkin->ensem.Vint = V_lj612;
  napkin->ensem.Vend = V_zero;
  napkin->ensem.Vext = Vext_harm;

  Rr_circlatt(napkin);
  Rend_close(napkin);
  // Rend_open(napkin);
  Rm_const(napkin, 1.0);

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

  for (size_t istep = 0;
      istep < napkin->ensem.nstep.thrm + napkin->ensem.nstep.prod;
      ++istep) {
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

    if (istep < napkin->ensem.nstep.thrm) {
      if (napkin->ensem.nstep.prodrec * napkin->ensem.istep.thrm >
          napkin->ensem.nstep.prod * napkin->ensem.istep.thrmrec) {
        disp_drift(napkin, driftfp);

        ++napkin->ensem.istep.thrmrec;
      }

      ++napkin->ensem.istep.thrm;
    } else {
      stats_accum(napkin->tde, est_tde(&napkin->ensem));

      if (napkin->ensem.nstep.prodrec * napkin->ensem.istep.prod >
          napkin->ensem.nstep.prod * napkin->ensem.istep.prodrec) {
        disp_drift(napkin, driftfp);
        disp_tde(napkin, tdefp);

        ++napkin->ensem.istep.prodrec;
      }

      ++napkin->ensem.istep.prod;
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
