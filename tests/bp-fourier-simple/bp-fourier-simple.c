/* A first attempt at basis pursuit */

#include <cm/sndfio.h>
#include <stdio.h>
#include <stdlib.h> 
#include <sndfile.h>
#include <cm_dict_gen.h>
#include <glpk.h> 
#include <assert.h> 
#include <signal.h>

/* To test the library we're making, this one uses a signal built from members
 * of the dictionary.
 * 
 * The goal is to solve the linear program:
 *
 * minimize |a|-1 (the l-1 norm o vector a)
 *
 * subject to
 * (phi)(a) = s
 *
 * where phi is a dictionary of waveforms
 * and s is the signal
 *
 */

int main(int argc, char **argv)
{

    /* Args are bp-fourier-proc dict num overcompleteness-parameter-l*/
    if (argc != 4) {
        fprintf(stderr, "args are %s which-dictionary-element "
                "overcompleteness-parameter-l signal-length\n", argv[0]);
        exit(1);
    }

    size_t k = atoi(argv[1]); /* which member of dictionary to use as signal */
    size_t l = atoi(argv[2]); /* How overcomplete the dictionary is */
    size_t n = atoi(argv[3]); /* How long the signal is */
    size_t p = n*l;           /* How many waveforms in the dictionary */
    size_t i, j;              /* row, column indices */




    /* Allocate fourier dictionary with enough waveforms */
    cm_dg_dict_t *phi = cm_dg_dict_alloc(n, p, CM_DG_DICT_FOURIER, NULL);
    cm_dg_fill_fourier(phi->data, n, l);

    /* copy a dictionary element into the signal for testing */
    double *s = (double*)malloc(sizeof(double)*n);
    memset(s, 0, sizeof(double)*n);
    for(i = 0; i < n; i++)
        s[i] = cm_dg_get_c(phi, i, k);

    /* Now we set up the linear program */

    /* Row indices of matrix phi */
    int *phi_i, *phi_j;
    /* entries of constraint matrix */
    double *A;
    assert((phi_i = (int*)malloc(sizeof(int)*(n*p + 1))));
    assert((phi_j = (int*)malloc(sizeof(int)*(n*p + 1))));
    assert((A = (double*)malloc(sizeof(double)*(n*p + 1))));

    /* create problem */
    glp_prob *lp = glp_create_prob();

    /* make into minimization problem */
    glp_set_obj_dir(lp, GLP_MIN);

    /* add rows for the signal vector */
    glp_add_rows(lp, n);
    /* add colums for constraints */
    glp_add_cols(lp, n);

    for (i = 0; i < n; i++) {
        /* nth row should equal the nth sample in the signal */
        glp_set_row_bnds(lp, i+1, GLP_FX, s[i], s[i]);
        /* the constraint on x is that it be greater or equal to 0 */
        glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
        /* each coefficient of the l-1 norm is just 1 */
        glp_set_obj_coef(lp, i+1, 1.);
    }

    /* Initialize constraint matrix */

    for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++){
            phi_i[i*p + j + 1] = i + 1;
            phi_j[i*p + j + 1] = j + 1;
            A[i*p + j + 1]   = cm_dg_get_c(phi, i, j);
        }
    }

    /* Load in the constraint matrix */
    glp_load_matrix(lp, n*p, phi_i, phi_j, A);

    /* Set up parameters for solving the problem */
    glp_smcp params;
    glp_init_smcp(&params);
    params.meth = GLP_DUALP;

    /* Solve the linear program */
    int err;
    if((err = glp_simplex(lp, &params))){
        fprintf(stderr, "Error %d in solving the LP.\n", err);
        exit(1);
    }

    for (i = 0; i < p; i++)
        printf("Column %d: %g\n", (int)i+1, glp_get_col_prim(lp,i+1));

    glp_delete_prob(lp);
    free(s);
    free(phi_i);
    free(phi_j);
    free(A);
    cm_dg_dict_free(phi);

    return(0);

}





