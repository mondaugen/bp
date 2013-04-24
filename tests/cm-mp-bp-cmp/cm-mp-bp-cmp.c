/* Test the loading of soundfiles into a dictionary */

#include <cm/sndfio.h> 
#include <cm_dict_gen.h>
#include <cm_bp.h> 
#include <cm_mp.h> 
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h> 
#include <assert.h> 

int
main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr, "Args are %s number-of-iterations normalize={1,0} "
                "signal-length LP-verbose={1,0}\n", argv[0]);
        exit(1);
    }

    size_t i,j,n,p,l;

    /* Length of signal */
    n = atoi(argv[3]);

    /* Using Fourier Dictionary, so l is how over-determined dictionary is */
    l = 2;

    /* Number of atoms in dictionary */
    p = n*l;

    /* make a fourier dictionary */
    cm_dg_dict_t *dict
        = cm_dg_dict_alloc(n, p, CM_DG_DICT_FOURIER, NULL);
    cm_dg_fill_fourier(dict->data, dict->n, l);

    /* Number of atoms in the signal construction */
    size_t gammaSize = 2;

    /* Columns of dictionary we use to construct the waveform */
    size_t Gamma[] =    {2,4};

    /* Scalars applying to the above columns */
    double Lambda[] =   {0.5, 0.75};


    printf("Fourier dictionary.\n");

    /* Print dictionary out */
    for (i = 0; i < dict->n; i++) {
        for (j = 0; j < dict->p; j++)
            printf("%+1.3f ", cm_dg_get_c(dict, i, j)); 
        printf("\n");
    }
    /* Indicate which atoms make up signal */
    printf("Coefficients of atoms that make up signal.\n");
    for (i = 0; i < dict->p; i++) {
        size_t k;
        int idx_in_gamma = 0;
        for (k = 0; k < gammaSize; k++) {
            if (Gamma[k] == i) {
                idx_in_gamma = 1;
                printf("%+1.3f ", Lambda[k]);
                break;
            }
        }
        if (!idx_in_gamma)
            printf("       ");
    }
    printf("\n\n");

    /* make a signal */
    double *s;
    assert((s = (double*)malloc(sizeof(double)*n)));
    double sigmax = 0;
    memset(s,0,sizeof(double)*n);
    for (i = 0; i < n; i++) {
        size_t k;
        for (k = 0; k < gammaSize; k++) {
            s[i] += Lambda[k] * cm_dg_get_c(dict, i, Gamma[k]);
            /* Keep track of signal maximum to normalize later */
            if (s[i] > sigmax)
                sigmax = s[j]; 
        }
    }

    if (*(argv[2]) == 1) {
        /* Normalize signal */
        for (i = 0; i < n; i++)
            s[i] /= sigmax;
        printf("Normalizing signal...\n");
    }
    else
        printf("Not normalizing signal...\n");

    printf("Signal to decompose\n");
    for (i = 0; i < n; i++)
        printf("%+1.3f ", s[i]);
    printf("\n\n");

    /* space for residual */
    double *r;
    assert((r = (double*)malloc(sizeof(double)*n)));

    /* space for residual + decomposed signal */
    double *mpcheck;
    assert((mpcheck = (double*)malloc(sizeof(double)*n)));

    /* coefficients go in here */
    double *x;
    assert((x = (double*)malloc(sizeof(double)*p)));

    /* Number of iterations */
    size_t K;

    K = atoi(argv[1]);

    printf("Doing %lu iterations of MP.\n", K);

    /* do MP on the signal with the constraint matrix dict */
    int err;
    if((err = cm_mp_decomp(x, r, dict->data, s, dict->n,
                    dict->p, K)))
        fprintf(stderr,"Error %d solving MP.\n", err);

    printf("Coefficient vector\n");
    for (i = 0; i < p; i++)
        printf("%+1.3f ", x[i]);
    printf("\n\n");

    printf("Residual\n");
    for (i = 0; i < n; i++)
        printf("%+1.3f ", r[i]);
    printf("\n\n");


    /* Check that matching pursuit was correct */
    memcpy(mpcheck,r,n*sizeof(double));
    /* vector view of mpcheck and coef vector x */
    gsl_vector_view mpc_view = gsl_vector_view_array(mpcheck,n);
    gsl_vector_view x_view = gsl_vector_view_array(x,p);
    /* matrix view of dictionary */
    gsl_matrix_view D_view = gsl_matrix_view_array(dict->data,dict->n,dict->p);
    /* calculate Dx + r, which should equal s */
    gsl_blas_dgemv(CblasNoTrans,1.,&D_view.matrix,&x_view.vector,
            1.,&mpc_view.vector);

    printf("Residual + Dx\n");
    for (i = 0; i < n; i++)
        printf("%+1.3f ", mpcheck[i]);
    printf("\n\n");

    printf("Doing LP.\n");

    /* do LP on the signal with constraint matrix dict */
    unsigned int lpargs = 0;
    if (atoi(argv[4]) == 1)
        lpargs |= CM_BP_LP_VERBOSE;
    if ((err = cm_bp_lp_solve(dict->data, dict->n, dict->p, x, s, lpargs)))
        fprintf(stderr,"Error %d solving LP.\n", err);

    printf("Coefficient vector\n");
    for (i = 0; i < p; i++)
        printf("%+1.3f ", x[i]);
    printf("\n\n");

    cm_dg_dict_free(dict);
    free(s);
    free(r);
    free(x);
    free(mpcheck);

    return(0);

}

