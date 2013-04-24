#ifndef CM_PDLBLP_H
#define CM_PDLBLP_H 

#include <stddef.h> 
#include <float.h> 
#include <limits.h> 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h> 
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_math.h> 

/* Multiply each row of matrix M by an element of vector v
 * If V is a matrix with v on the main diagonal, then the result is M = VM */
int
_matrix_rows_mul(gsl_matrix *M, gsl_vector *v);

/* find max{ p : x + p*dx >= 0 } */
int
_p_max_find(double *p, gsl_vector *x, gsl_vector *dx);

/* solve min{c^Tx + 1/2*|gamma*x|^2 + 1/2|p|^2} subject to Ax + delta*p = b,
 * x >= 0.
 * For BP, c is usually e (a vector of all ones), A(n,p) is the dictionary of p
 * waveforms of length n, b(n) is the signal and the resulting x(p) is the vector of
 * coefficients characterizing the basis. Gamma and delta are the regularization
 * parameters "feasibility tolerance" and "duality gap tolerance", K is the
 * maximum number of iterations of the algorithm to carry out.
 */
int
pdlblp_solve(gsl_vector *x, gsl_matrix *A, gsl_vector *b, gsl_vector *c,
        double gamma, double delta, size_t K, double yinit,
	double zinit, double muinit);

#endif /* CM_PDLBLP_H */
