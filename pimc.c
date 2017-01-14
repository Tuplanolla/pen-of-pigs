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

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  char const* const shortstr = "s:d:N:M:k:K:h:p:H:P:L:m:T:";
  char const* const longstr[] = {
    "sys", "ndim", "npoly", "nbead", "nsubdiv", "ndiv",
    "nthrm", "nprod", "nthrmrec", "nprodrec",
    "length", "mass", "temp"
  };
  char const* const sys[] = {
    "qho"
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
  double temp = 1.0;

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
      case 'T':
        if (opt_parse_fp(&temp, 0.0, INFINITY))
          continue;
        break;
    }

    (void) fprintf(stderr, "Failed to parse argument list.\n");

    return EXIT_FAILURE;
  }

  switch (isys) {
    case 0:
      {
        double const beta = 1.0 / temp;

        struct sim* const sim = sim_alloc(ndim, npoly, nbead, nsubdiv, ndiv,
            nthrm, nprod, nthrmrec, nprodrec,
            true, false, length, 1.0, beta);
        if (sim == NULL) {
          (void) fprintf(stderr, "Failed to allocate memory.\n");

          return EXIT_FAILURE;
        }

        double const q = omega / 2.0;
        double const e = (double) ndim * q / tanh(q * beta);
        (void) printf("Expected for QHO: E = %f (T = %f)\n", e, temp);

        sim_set_potext(sim, potext_harm);

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
