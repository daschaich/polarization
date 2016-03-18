// -----------------------------------------------------------------
// Defines and subroutine declarations SU(N) polarization studies
#ifndef _SUN_H
#define _SUN_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge group is SU(NCOL)
// Assume only fundamental fermions with NCOL = NCOL
#define NCOL 2

#if (NCOL != 3)
#ifdef FAST
  #error "FAST only works if NCOL=3!"
#endif
#endif

typedef struct { fcomplex e[NCOL][NCOL]; } fsu3_matrix;
typedef struct { fcomplex c[NCOL]; } fsu3_vector;

// Anti-hermitian matrices for general NCOL
typedef struct {
  fcomplex m[NCOL * (NCOL - 1) / 2];
  float im_diag[NCOL];
} fanti_hermitmat;

typedef struct { dcomplex e[NCOL][NCOL]; } dsu3_matrix;
typedef struct { dcomplex c[NCOL]; } dsu3_vector;
typedef struct {
  dcomplex m[NCOL * (NCOL - 1) / 2];
  double im_diag[NCOL];
} danti_hermitmat;

#if (PRECISION==1)
#define su3_matrix      fsu3_matrix
#define su3_vector      fsu3_vector
#define anti_hermitmat  fanti_hermitmat
#else
#define su3_matrix      dsu3_matrix
#define su3_vector      dsu3_vector
#define anti_hermitmat  danti_hermitmat
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Matrix operations
// In file trace_su3.c
//complex trace_su3(su3_matrix *a);

// In file realtr.c
Real realtrace_su3(su3_matrix *a, su3_matrix *b);

// In file complextr.c
//complex complextrace_su3(su3_matrix *a, su3_matrix *b);

// In file addmat.c
//void add_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file submat.c
//void sub_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file s_m_mat.c
//void scalar_mult_su3_matrix(su3_matrix *src, Real scalar, su3_matrix *dest);

// In file s_m_a_mat.c
//void scalar_mult_add_su3_matrix(su3_matrix *src1, su3_matrix *src2,
//                                Real scalar, su3_matrix *dest);

// In file s_m_s_mat.c
//void scalar_mult_sub_su3_matrix(su3_matrix *src1, su3_matrix *src2,
//                                Real scalar, su3_matrix *dest);

// In file cs_m_mat.c
//void c_scalar_mult_su3mat(su3_matrix *src, complex *scalar,
//                          su3_matrix *dest);

// In file cs_m_a_mat.c
//void c_scalar_mult_add_su3mat(su3_matrix *src1, su3_matrix *src2,
//                              complex *scalar, su3_matrix *dest);

// In file cs_m_s_mat.c
//void c_scalar_mult_sub_su3mat(su3_matrix *src1, su3_matrix *src2,
//                              complex *scalar, su3_matrix *dest);

// In file su3_adjoint.c
//void su3_adjoint(su3_matrix *a, su3_matrix *b);

// In file clear_mat.c
//void clear_su3mat(su3_matrix *dest);

// In file su3mat_copy.c
//void su3mat_copy(su3_matrix *a, su3_matrix *b);

// In file dumpmat.c
//void dumpmat(su3_matrix *m);

// In file m_mat_nn.c
void mult_su3_nn(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file m_mat_na.c
//void mult_su3_na(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file m_mat_an.c
void mult_su3_an(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// c <-- a * b, in file m_matvec.c
//void mult_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- c + a * b, in file m_matvec_s.c
//void mult_su3_mat_vec_sum(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- c - a * b, in file m_matvec_ns.c
//void mult_su3_mat_vec_nsum(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- adag * b, in file m_amatvec.c
//void mult_adj_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- c + adag * b, in file m_amatvec_s.c
//void mult_adj_su3_mat_vec_sum(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- c - adag * b, in file m_amatvec_ns.c
//void mult_adj_su3_mat_vec_nsum(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- b * a, in file m_vecmat.c
//void mult_su3_vec_mat(su3_vector *b, su3_matrix *a, su3_vector *c);

// c <-- b * adag, in file m_vecamat.c
//void mult_su3_vec_adj_mat(su3_vector *b, su3_matrix *a, su3_vector *c);

// c <-- c + b * adag, in file m_vecamat_s.c
//void mult_su3_vec_adj_mat_sum(su3_vector *b, su3_matrix *a, su3_vector *c);

// In file m_amv_4dir.c
//void mult_adj_su3_mat_vec_4dir(su3_matrix *a, su3_vector *b, su3_vector *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Miscellaneous routines
// In file gaussrand.c
Real gaussian_rand_no(double_prn *prn_pt);

#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

#endif
// -----------------------------------------------------------------
