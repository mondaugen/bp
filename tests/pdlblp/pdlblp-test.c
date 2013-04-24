/* Test of primal-dual log-barrier LP algorithm */

#include <cm_pdlblp.h>
#include <stddef.h> 
#include <stdio.h> 

int
main(void)
{
    /* Trivial test */

    /* Regularization parameters */
    double gamma = 1.E-4;
    double delta = 0.;

    /* Number of iterations */
    size_t K = 2;

    /* Value to initialize to */
    double init = 0.1;

    /* Dictionary/constraint matrix */
    gsl_matrix *A = gsl_matrix_alloc(4,4);
    gsl_matrix_set_identity(A);

    /* Coefficient vector */
    gsl_vector *x = gsl_vector_alloc(4);
    /* Give a feasible solution */
    gsl_vector_set(x,0,1.);
    gsl_vector_set(x,1,2.);
    gsl_vector_set(x,2,3.);
    gsl_vector_set(x,3,4.);

    /* Signal vector */
    gsl_vector *b = gsl_vector_alloc(4);
    gsl_vector_set(b,0,1.);
    gsl_vector_set(b,1,2.);
    gsl_vector_set(b,2,3.);
    gsl_vector_set(b,3,4.);

    /* Coefficient vector scalars */
    gsl_vector *c = gsl_vector_alloc(4);
    gsl_vector_set_all(c,1.);

    pdlblp_solve(x, A, b, c, 1.E-4, 1.E-4, K, 0., init, 0.001);

    gsl_vector_fprintf(stdout, x, "%f");

    gsl_matrix_free(A);
    gsl_vector_free(x);
    gsl_vector_free(b);
    gsl_vector_free(c);

    return(0);

}

