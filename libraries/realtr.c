// -----------------------------------------------------------------
// Return real trace of matrix product adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real realtrace_su3(su3_matrix *a, su3_matrix *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      sum += a->e[i][j].real * b->e[i][j].real
           + a->e[i][j].imag * b->e[i][j].imag;
    }
  }
  return sum;
}
// -----------------------------------------------------------------
