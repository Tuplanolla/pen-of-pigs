#include "exts.h"
#include "lims.h"
#include "opt.h"
#include "sim.h"
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double const omega = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double potext_harm(struct ens const* const ens,
    struct bead const* const r) {
  return gsl_pow_2(omega) * ens_norm2(ens, r) / 2.0;
}

static double const epsilon = 10.22;
static double const sigma = 2.556;

__attribute__ ((__nonnull__, __pure__))
static double pot_lj612(struct ens const* const ens,
    struct bead const* const r0, struct bead const* const r1) {
  double const d = ens_dist2(ens, r0, r1);

  if (d == 0.0)
    return INFINITY;
  else {
    double const sigmad2 = gsl_pow_2(sigma) / d;

    return 4.0 * epsilon * (gsl_pow_6(sigmad2) - gsl_pow_3(sigmad2));
  }
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  char const* const shortstr = "s:d:N:M:k:K:h:p:H:P:L:m:t:";
  char const* const longstr[] = {
    "sys", "ndim", "npoly", "nbead", "nsubdiv", "ndiv",
    "nthrm", "nprod", "nthrmrec", "nprodrec",
    "mass", "length", "imagtime"
  };
  char const* const sys[] = {
    "qho", "he4"
  };

  size_t isys = 0;
  size_t ndim = 1;
  size_t npoly = 1;
  size_t nbead = 1;
  size_t nsubdiv = 1;
  size_t ndiv = 1;
  size_t nthrm = 0;
  size_t nprod = 0;
  size_t nthrmrec = 0;
  size_t nprodrec = 0;
  double length = 1.0;
  double mass = 1.0;
  double tau = 1.0;

  for ever {
    int const i = opt_parse(argc, argv, shortstr, longstr);
    if (i == -1)
      break;

    switch ((char) i) {
      case 's':
        if (opt_parse_str(&isys, sys, sizeof sys / sizeof *sys))
          continue;
        break;
      case 'd':
        if (opt_parse_size(&ndim, 1, DIM_MAX))
          continue;
        break;
      case 'N':
        if (opt_parse_size(&npoly, 1, POLY_MAX))
          continue;
        break;
      case 'M':
        if (opt_parse_size(&nbead, 1, BEAD_MAX))
          continue;
        break;
      case 'k':
        if (opt_parse_size(&nsubdiv, 1, SUBDIV_MAX))
          continue;
        break;
      case 'K':
        if (opt_parse_size(&ndiv, 1, DIV_MAX))
          continue;
        break;
      case 'h':
        if (opt_parse_size(&nthrm, 0, STEP_MAX))
          continue;
        break;
      case 'p':
        if (opt_parse_size(&nprod, 0, STEP_MAX))
          continue;
        break;
      case 'H':
        if (opt_parse_size(&nthrmrec, 0, STEP_MAX))
          continue;
        break;
      case 'P':
        if (opt_parse_size(&nprodrec, 0, STEP_MAX))
          continue;
        break;
      case 'L':
        if (opt_parse_fp(&length, 0.0, INFINITY))
          continue;
        break;
      case 'm':
        if (opt_parse_fp(&mass, 0.0, INFINITY))
          continue;
        break;
      case 't':
        if (opt_parse_fp(&tau, 0.0, INFINITY))
          continue;
        break;
    }

    (void) fprintf(stderr, "Failed to parse argument list.\n");

    return EXIT_FAILURE;
  }

  switch (isys) {
    case 0:
      {
        double const q = omega / 2.0;
        double const e = (double) ndim * q;
        (void) printf("Expected for QHO: E = %f\n", e);

        struct sim* const sim = sim_alloc(ndim, npoly, nbead, nsubdiv, ndiv,
            nthrm, nprod, nthrmrec, nprodrec,
            false, length, mass, tau);
        if (sim == NULL) {
          (void) fprintf(stderr, "Failed to allocate memory.\n");

          return EXIT_FAILURE;
        }

        sim_perm_open(sim, NULL);
        sim_set_potext(sim, potext_harm);

        if (!sim_run(sim)) {
          (void) fprintf(stderr, "Failed to run simulation.\n");

          return EXIT_FAILURE;
        }

        sim_free(sim);

        return EXIT_SUCCESS;
      }
    case 1:
      {
        // Equilibrium conditions:
        // d = 3
        // N = 64 [16]
        // M = 15
        // epsilon = 0.0125 K^-1 (tau = 0.2 K^-1)
        // rho = 0.02186 A^-3 (L = 14.306 A [9.0 A])
        double const en = -7.32;
        (void) printf("Expected for He-4: E = %f (E / N = %f)\n", en * npoly, en);

        struct sim* const sim = sim_alloc(ndim, npoly, nbead, nsubdiv, ndiv,
            nthrm, nprod, nthrmrec, nprodrec,
            true, length, 1.0, tau);
        if (sim == NULL) {
          (void) fprintf(stderr, "Failed to allocate memory.\n");

          return EXIT_FAILURE;
        }

        sim_perm_open(sim, NULL);
        sim_set_potint(sim, pot_lj612);
        sim_place_lattice(sim, sim_placer_random, NULL);

        if (!sim_run(sim)) {
          (void) fprintf(stderr, "Failed to run simulation.\n");

          return EXIT_FAILURE;
        }

        sim_free(sim);

        return EXIT_SUCCESS;
      }
    default:
      return EXIT_FAILURE;
  }
}
