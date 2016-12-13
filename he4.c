#include "exts.h"
#include "lims.h"
#include "opt.h"
#include "sim.h"
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double const epsilon = 2.0;
static double const sigma = 0.2;

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
  char const* const shortstr = "d:N:M:K:h:p:H:P:L:";
  char const* const longstr[] = {
    "ndim", "npoly", "nbead", "nsubdiv",
    "nthrm", "nprod", "nthrmrec", "nprodrec",
    "length"
  };

  size_t ndim = 1;
  size_t npoly = 1;
  size_t nbead = 1;
  size_t nsubdiv = 1;
  size_t nthrm = 0;
  size_t nprod = 0;
  size_t nthrmrec = 0;
  size_t nprodrec = 0;
  double L = 1.0;

  for ever {
    int const i = opt_parse(argc, argv, shortstr, longstr);
    if (i == -1)
      break;

    switch ((char) i) {
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
      case 'K':
        if (opt_parse_size(&nsubdiv, 1, SUBDIV_MAX))
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
        if (opt_parse_fp(&L, 0.0, INFINITY))
          continue;
        break;
    }

    (void) fprintf(stderr, "Failed to parse argument list.\n");

    return EXIT_FAILURE;
  }

  struct sim* const sim = sim_alloc(ndim, npoly, nbead, nsubdiv,
      nthrm, nprod, nthrmrec, nprodrec,
      true, L, 1.0, 100.0);
  if (sim == NULL) {
    (void) fprintf(stderr, "Failed to allocate memory.\n");

    return EXIT_FAILURE;
  }

  sim_set_potint(sim, pot_lj612);
  sim_perm_open(sim, NULL);
  sim_place_random(sim, sim_placer_point, NULL);

  if (!sim_run(sim)) {
    (void) fprintf(stderr, "Failed to run simulation.\n");

    return EXIT_FAILURE;
  }

  sim_free(sim);

  return EXIT_SUCCESS;
}
