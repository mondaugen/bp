/* A first attempt at basis pursuit */

#include <cm/sndfio.h>
#include <stdio.h>
#include <stdlib.h> 
#include <sndfile.h>
#include <cm_dict_gen.h>
#include <glpk.h> 
#include <assert.h> 
#include <signal.h>
#include <math.h> 

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
    if (argc != 7) {
        fprintf(stderr, "args are %s overcompleteness-parameter-l"
               " signal-length frequency amplitude phase sin|cos\n", argv[0]);
        exit(1);
    }

    size_t l = atoi(argv[1]); /* How overcomplete the dictionary is */
    size_t n = atoi(argv[2]); /* How long the signal is */
    size_t p = n*l;           /* How many waveforms in the dictionary */
    size_t i, j;              /* row, column indices */
    double f = atof(argv[3]); /* frequency of signal */
    double a = atof(argv[4]); /* amplitude of signal */
    double ph = atof(argv[5]); /* phase of signal */
    int    v = atoi(argv[6]); /* if non-zero, use sine, otherwise cosine */




    /* Allocate fourier dictionary with enough waveforms */
    cm_dg_dict_t *phi = cm_dg_dict_alloc(n, p, CM_DG_DICT_FOURIER, NULL);
    cm_dg_fill_fourier(phi->data, n, l);

    /* Prepare a signal for testing.
     * The signal is a sinusoid
     * frequency is normalized by signal length, so f=1 fits one period of wave
     * in the signal, f=2 fits 2 etc.
     * phase ph is 0-1 (phase offset of wave)
     * a is the amplitude
     */
    double *s = (double*)malloc(sizeof(double)*n);
    memset(s, 0, sizeof(double)*n);
    for (i = 0; i < n; i++) {
        if( v )
            s[i] = a * sin(2. * M_PI * (f * i / (double)n + ph));
        else
            s[i] = a * cos(2. * M_PI * (f * i / (double)n + ph));
    }

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
        /* x is a free variable */
        glp_set_col_bnds(lp, i+1, GLP_FR, 0.0, 0.0);
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





