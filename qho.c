#include "sim.h"
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__attribute__ ((__const__, __nonnull__))
static double V_zero(
    __attribute__ ((__unused__)) struct ensem const* const ensem,
    __attribute__ ((__unused__)) struct bead const* const r0,
    __attribute__ ((__unused__)) struct bead const* const r1) {
  return 0.0;
}

__attribute__ ((__nonnull__, __const__))
static double Vext_zero(
    __attribute__ ((__unused__)) struct ensem const* const ensem,
    __attribute__ ((__unused__)) struct bead const* const r) {
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

static double const omega = 1.0;

__attribute__ ((__nonnull__, __pure__))
static double Vext_harm(struct ensem const* const ensem,
    struct bead const* const r) {
  return gsl_pow_2(omega) * bead_norm2(ensem, r) / 2.0;
}

int main(void) {
  size_t const ndim = 3;
  size_t const npoly = 3;
  size_t const nbead = 8;
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

  sim_run(ndim, npoly, nbead, nsubdiv, nthrm, nprod, nthrmrec, nprodrec,
      true, 10.0, beta, V_lj612, V_zero, Vext_harm);

  return EXIT_SUCCESS;
}
