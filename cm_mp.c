/* Matching pursuit algorithms */

#include <cm_mp.h>
#include <string.h> 
#include <math.h> 

/* Do matching pursuit on signal s with atoms in D. Do K iterations. Result goes
 * into x, residual in r.
 * s, r length n
 * x length p
 * D has n rows and p columns
 * Note: x gets initialized to 0.
 */
int
cm_mp_decomp(double *x, double *r, double *D, double *s, size_t n,
        size_t p, size_t K)
{
    
    /* Initialize x to 0 */
    memset(x, 0, sizeof(double)*p);

    /* Initialize residual as signal */
    memcpy(r, s, sizeof(double)*n);

    /* Make vector-view of r and matrix-view of D */
    gsl_vector_view r_view = gsl_vector_view_array(r,n);
    gsl_matrix_view D_view = gsl_matrix_view_array(D,n,p);

    size_t i, j;

    while (K-- > 0) {

        /* The current vector and the vector with the greatest inner-product*/
        gsl_vector_view c, gv;

        /* The greatest inner-product */
        double gip = 0;

        /* The index of the greatest inner-product */
        size_t gi;

        for (i = 0; i < p; i++) {
            c = gsl_matrix_column(&D_view.matrix, i);
            double ip;
            gsl_blas_ddot(&r_view.vector, &c.vector, &ip);
            /* If the absolute value of this inner-product greater than the
             * absolute value of the last greatest, store this inner-product */
            if ( fabs(ip) > fabs(gip) ) {
                gip = ip;
                gi = i;
                gv = c;
            }
        }

        if (gip != 0.) {
            
            /* Add contrubution of the inner-product to the corresponding
             * coordinate */
            x[gi] += gip;

            /* Subtract contribution of this vector from the residual */
            gsl_blas_daxpy(-1. * gip, &gv.vector, &r_view.vector);

        }

    }

    return(0);

}









