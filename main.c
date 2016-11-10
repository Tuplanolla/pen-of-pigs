#include "he4.h"
#include "qho.h"
#include "report.h"
#include <stdlib.h>

#include "floating.h"
#include <stdio.h>

/* int main(int const n, char** const xs) { */
int main(void) {
  reset();
  warn(reset);

  bool const p =
#ifdef PEN_OF_HE4
    he4() &&
#endif
#ifdef PEN_OF_QHO
    qho() &&
#endif
    true;

  /*
  for (double i = -8; i <= 8; i += 1)
    printf("%f\n", wrap(i, 3));
  */

  return p ? EXIT_SUCCESS : EXIT_FAILURE;
}
