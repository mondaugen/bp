/* A first attempt at basis pursuit */

#include <cm/sndfio.h>
#include <stdio.h>
#include <stdlib.h> 
#include <sndfile.h>
#include <cm_dict_gen.h>
#include <glpk.h> 
#include <assert.h> 

/* The goal is to solve the linear program:
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

    /* Args are bp-fourier-proc sndfile overcompleteness-parameter-l*/
    if (argc != 6) {
        fprintf(stderr, "args are %s path-to-soundfile "
                "overcompleteness-parameter-l sampling-rate simplex|interior "
		"dictionary\n", argv[0]);
        exit(1);
    }

    SF_INFO sfinfo;

    double *s;
    
    if(!(s = cm_load_sf_to_dbl_ary(argv[1], &sfinfo))) {
        fprintf(stderr,"Error reading soundfile. %s\n", argv[1]);
        exit(1);
    }

    if (sfinfo.channels != 1) {
        fprintf(stderr,"Mono files only please.\n");
        free(s);
        exit(1);
    }

    size_t l = atoi(argv[2]); /* How overcomplete the dictionary is */
    size_t n = sfinfo.frames; /* How long the signal is */
    size_t p = n*l;           /* How many waveforms in the dictionary */
    double sr= atof(argv[3]);
    int simplex = atoi(argv[4]); /* if non-zero, simplex, else interior */
    char dict = *(argv[5]);

    /* Allocate dictionary */
    cm_dg_dict_t *phi;

    /* Now we set up the linear program */

    /* Row indices of matrix phi */
    int *phi_i, *phi_j;
    /* entries of constraint matrix */
    double *A;

    /* create problem */
    glp_prob *lp = glp_create_prob();

    /* make into minimization problem */
    glp_set_obj_dir(lp, GLP_MIN);

    size_t i, j;

    if (dict == 'f') {
	fprintf(stderr,"Using fourier dictionary.\n");

	/* use fourier dictionary */
	phi = cm_dg_dict_alloc(n, p, CM_DG_DICT_FOURIER, NULL);
	cm_dg_fill_fourier(phi->data, n, l);

	assert((phi_i = (int*)malloc(sizeof(int)*(n*p + 1))));
	assert((phi_j = (int*)malloc(sizeof(int)*(n*p + 1))));
	assert((A = (double*)malloc(sizeof(double)*(n*p + 1))));
	/* add rows for the signal vector */
	glp_add_rows(lp, n);
	/* add colums for constraints */
	glp_add_cols(lp, p);

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

    } else if (dict == 'i') {
	fprintf(stderr,"Using identity dictionary.\n");

	phi = cm_dg_dict_alloc(n, n, CM_DG_DICT_IDENTITY, NULL);
	cm_dg_fill_identity(phi->data, n);

	/* Use identity matrix as dictionary */
	assert((phi_i = (int*)malloc(sizeof(int)*(n*n + 1))));
	assert((phi_j = (int*)malloc(sizeof(int)*(n*n + 1))));
	assert((A = (double*)malloc(sizeof(double)*(n*n + 1))));

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
	    for (j = 0; j < n; j++){
		phi_i[i*n + j + 1] = i + 1;
		phi_j[i*n + j + 1] = j + 1;
		A[i*n + j + 1]   = cm_dg_get_c(phi, i, j);
	    }
	}
	/* Load in the constraint matrix */
	glp_load_matrix(lp, n*n, phi_i, phi_j, A);

    } else {
	glp_delete_prob(lp);
	free(s);
	fprintf(stderr,"Dictionary argument accepts i|f\n");
	exit(1);
    }

    if (simplex) {

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

    } else {

        /* Set up parameters for solving the problem */
        glp_smcp params;
        glp_init_smcp(&params);
        params.meth = GLP_DUALP;

        /* Solve the linear program */
        int err;
        if((err = glp_interior(lp, NULL))){
            fprintf(stderr, "Error %d in solving the LP.\n", err);
            exit(1);
        }

    }

    if (dict == 'f') {

	for (i = 0; i <= (p/2); i++) {

	    printf("Cos frequency: %g: amplitude: %g\n", (double)i/(double)p*sr,
		    glp_get_col_prim(lp,i+1));

	}
	for ( ; i < p; i++) {

	    printf("Sin frequency: %g: amplitude: %g\n",
		    (double)(i-(p/2))/(double)p*sr,
		    glp_get_col_prim(lp,i+1));

	}
    } else if (dict == 'i') {

	for (i = 0; i < n; i++) {

	    printf("Vector %lu: amplitude: %g\n", i,
		    glp_get_col_prim(lp,i+1));

	}
    }

    glp_delete_prob(lp);
    free(s);
    free(phi_i);
    free(phi_j);
    free(A);
    cm_dg_dict_free(phi);

    return(0);

}





