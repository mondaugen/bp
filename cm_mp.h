#ifndef CM_MP_H
#define CM_MP_H 

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h> 
#include <stddef.h> 

/* Do matching pursuit on signal s with atoms in D. Do K iterations. Result goes
 * into x, residual in r.
 * s, r length n
 * x length p
 * D has n rows and p columns
 * Note: x gets initialized to 0.
 */
int
cm_mp_decomp(double *x, double *r, double *D, double *s, size_t n,
        size_t p, size_t K);

#endif /* CM_MP_H */
