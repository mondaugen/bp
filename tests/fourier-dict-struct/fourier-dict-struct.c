/* Testing the generation of fourier dictionaries */

#include <stddef.h> 
#include <stdlib.h> 
#include <stdio.h>
#include "cm_dict_gen.h" 

int
main(int argc, char **argv)
{
    /* args are l and n */
    if (argc != 3) {
        fprintf(stderr, "Args are %s l and n.\n", argv[0]);
        exit(1);
    }

    size_t l = atoi(argv[1]);
    size_t n = atoi(argv[2]);

    cm_dg_dict_t *fdict = cm_dg_dict_alloc(n, l*n, CM_DG_DICT_FOURIER, NULL);

    int err;
    if((err = cm_dg_fill_fourier(fdict->data, n , l))){
        fprintf(stderr,"Error %d filling with dictionary.\n", err);
        exit(1);
    }

    size_t i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < (l*n); j++)
            printf("%+1.1f ", cm_dg_get_c(fdict, i, j));
        printf("\n");
    }

    cm_dg_dict_free(fdict);

    exit(0);

}
