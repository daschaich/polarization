// -----------------------------------------------------------------
// Include files for SU(N) pure-gauge polarization studies
#include "../include/config.h"  // Keep this first
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>             // For print_var.c, setup.c, gauge_info.c
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Basic observables
void plaquette(double *ss_plaq, double *st_plaq);

// Other routines in library_util.c that loop over all sites
void gauge_field_copy_f(field_offset src, field_offset dest);
void shiftmat(su3_matrix_f *dat, su3_matrix_f *temp, int dir);
// -----------------------------------------------------------------
