#include "exts.h"
#include "sim.h"
#include <gsl/gsl_math.h>
#include <stdio.h>
#include <stdlib.h>

static double const epsilon = 4.0;
static double const sigma = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double pot_lj612(struct ensem const* const ensem,
    struct bead const* const r0, struct bead const* const r1) {
  double const sigmar2 = gsl_pow_2(sigma) / bead_dist2(ensem, r0, r1);

  return 4.0 * epsilon * (gsl_pow_6(sigmar2) - gsl_pow_3(sigmar2));
}

int main(int const n, char** const xs) {
  size_t const ndim = 2;
  size_t const npoly = 8;
  size_t const nbead = 8;
  size_t const nsubdiv = ndim == 1 ? 256 : ndim == 2 ? 16 : 4;
  size_t const nthrm = 1 << 14;
  size_t const nprod = 1 << 18;
  size_t const nthrmrec = 1 << 4;
  size_t const nprodrec = 1 << 8;

  double const T = 0.1;
  double const beta = 1.0 / T;

  return sim_run(ndim, npoly, nbead, nsubdiv, nthrm, nprod, nthrmrec, nprodrec,
      true, 10.0, beta, pot_lj612, pot_zero, potext_zero) ?
    EXIT_SUCCESS : EXIT_FAILURE;
}
