/* Test the loading of soundfiles into a dictionary */

#include <cm/sndfio.h> 
#include <cm_dict_gen.h>
#include <cm_bp.h> 
#include <cm_mp.h> 
#include <stdlib.h> 
#include <stdio.h>

int
main(int argc, char **argv)
{
    size_t i,j;

    /* make an identity dictionary */
    cm_dg_dict_t *id_dict
        = cm_dg_dict_alloc(5, 5, CM_DG_DICT_IDENTITY, NULL);
    cm_dg_fill_identity(id_dict->data, id_dict->n);

    printf("Identity dictionary.\n");
    /* Print dictionary out */
    for (i = 0; i < id_dict->n; i++) {
        for (j = 0; j < id_dict->p; j++)
            printf("%1.3f ", cm_dg_get_c(id_dict, i, j)); 
        printf("\n");
    }
    printf("\n");

    /* make a signal */
    double s[] = {5., 2., 3., 4., 1.};

    /* space for residual */
    double r[5];

    /* coefficients go in here */
    double x[5];

    /* Number of iterations */
    size_t K;

    if (argc == 2)
        K = atoi(argv[1]);
    else
        K = 5;

    printf("Doing %lu iterations.\n", K);

    /* do MP on the signal with the constraint matrix id_dict */
    int err;
    if((err = cm_mp_decomp(x, r, id_dict->data, s, id_dict->n,
                    id_dict->p, K)))
        fprintf(stderr,"Error %d solving MP.\n", err);

    printf("Coefficient vector\n");
    for (i = 0; i < 5; i++)
        printf("%1.3f ", x[i]);
    printf("\n");

    cm_dg_dict_free(id_dict);

    return(0);

}

