#include "he4.h"
#include "qho.h"
#include "report.h"
#include <stdlib.h>

/* int main(int const n, char** const xs) { */
int main(void) {
  reset();

  bool const p =
#ifdef PEN_OF_HE4
    he4() &&
#endif
#ifdef PEN_OF_QHO
    qho() &&
#endif
    true;

  return p ? EXIT_SUCCESS : EXIT_FAILURE;
}
