// -----------------------------------------------------------------
// Matrix multiplication with no adjoints
// c <-- a * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_su3_nn(su3_matrix *a, su3_matrix *b, su3_matrix *c) {
  register int i, j, k;
  register complex x, y;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      x.real = 0.0;
      x.imag = 0.0;
      for (k = 0; k < NCOL; k++) {
        CMUL(a->e[i][k], b->e[k][j], y);
        CSUM(x, y);
      }
      c->e[i][j] = x;
    }
  }
}
#else   // FAST version for NCOL=3 only
void mult_su3_nn(su3_matrix *a, su3_matrix *b, su3_matrix *c) {
  int i, j;
  register Real t, ar, ai, br, bi, cr, ci;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      ar = a->e[i][0].real;
      ai = a->e[i][0].imag;
      br = b->e[0][j].real;
      bi = b->e[0][j].imag;
      cr = ar * br;
      t = ai * bi;
      cr -= t;
      ci = ar * bi;
      t = ai * br;
      ci += t;

      ar = a->e[i][1].real;
      ai = a->e[i][1].imag;
      br = b->e[1][j].real;
      bi = b->e[1][j].imag;
      t = ar * br;
      cr += t;
      t = ai * bi;
      cr -= t;
      t = ar * bi;
      ci += t;
      t = ai * br;
      ci += t;

      ar = a->e[i][2].real;
      ai = a->e[i][2].imag;
      br = b->e[2][j].real;
      bi = b->e[2][j].imag;
      t = ar * br;
      cr += t;
      t = ai * bi;
      cr -= t;
      t = ar * bi;
      ci += t;
      t = ai * br;
      ci += t;

      c->e[i][j].real = cr;
      c->e[i][j].imag = ci;
    }
  }
}
#endif
// -----------------------------------------------------------------
