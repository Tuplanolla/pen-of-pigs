#include "floating.h"
#include "report.h"
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <unistd.h>

int main(int const off_by_one, char** const mutable_storage) {
  reset();
  warn(main);
  unsigned int ok = sleep(1);
  warn(main);

  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
