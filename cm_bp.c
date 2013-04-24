#include <stdlib.h>
#include <stdio.h> 
#include <assert.h> 
#include <glpk.h> 

/* Basis pursuit algorithms */

/* Given a dictionary D, of n rows and p columns (n is signal length, p is
 * number of items in the dictionary), a signal s of length n and a coefficient
 * vector x of length p, solve the linear program min|x| subject to Dx = s.
 * D is stored in row-major order. */
int
cm_bp_lp_solve(double *D, size_t n, size_t p, double *x, double *s,
        int solveflags)
{
    int err = 0;
    size_t i, j;

    /* create linear programming problem */
    glp_prob *lp = glp_create_prob();

    /* make into minimization problem */
    glp_set_obj_dir(lp, GLP_MIN);

    /* Linear program must have its own constraint matrix and index matrices, so
     * allocate these */

    /* Row indices of matrix phi */
    int *phi_i, *phi_j;

    /* entries of constraint matrix */
    double *A;
	assert((phi_i = (int*)malloc(sizeof(int)*(n*p + 1))));
	assert((phi_j = (int*)malloc(sizeof(int)*(n*p + 1))));
	assert((A = (double*)malloc(sizeof(double)*(n*p + 1))));

	/* add rows for the signal vector */
	glp_add_rows(lp, n);
	
    /* add colums for constraints */
	glp_add_cols(lp, p);

    /* set up specific conditions for BP */
	for (i = 0; i < n; i++) {
	    /* nth row should equal the nth sample in the signal */
	    glp_set_row_bnds(lp, i+1, GLP_FX, s[i], s[i]);
	    /* x is a free variable */
	    glp_set_col_bnds(lp, i+1, GLP_FR, 0.0, 0.0);
	    /* each coefficient of the l-1 norm is just 1 */
	    glp_set_obj_coef(lp, i+1, 1.);
	}

	/* Initialize constraint matrix (indexed starting at 1, unlike c)*/
	for (i = 0; i < n; i++) {
	    for (j = 0; j < p; j++){
		phi_i[i*p + j + 1] = i + 1;
		phi_j[i*p + j + 1] = j + 1;
		A[i*p + j + 1]   = D[i*p + j];
	    }
	}

	/* Load in the constraint matrix */
	glp_load_matrix(lp, n*p, phi_i, phi_j, A);

    /* Only solves with simplex algorithm */

    /* Set up parameters for solving the problem */
    glp_smcp params;
    glp_init_smcp(&params);
    params.meth = GLP_DUALP;

    /* Solve the linear program */
    if((err = glp_simplex(lp, &params)) == 0) {
        /* Put result into x */
        for (i = 0; i < p; i++)
            x[i] = glp_get_col_prim(lp, i+1);
    }

    glp_delete_prob(lp);
    free(phi_i);
    free(phi_j);
    free(A);

    return(err);

}








