/* Test the loading of soundfiles into a dictionary */

#include <cm/sndfio.h> 
#include <cm_dict_gen.h>
#include <stdlib.h> 
#include <stdio.h>

int
main(void)
{
    size_t i,j;

    /* Load the junk from soundfiles */
    cm_dg_dict_t *dict = cm_dg_dict_alloc_from_dir(10,"./","*.wav");

    printf("Soundfile dictionary.\n");
    /* Print dictionary out */
    for (i = 0; i < dict->n; i++) {
        for (j = 0; j < dict->p; j++)
            printf("%1.3f ", cm_dg_get_c(dict, i, j)); 
        printf("\n");
    }

    /* Also make an identity dictionary */
    cm_dg_dict_t *id_dict
        = cm_dg_dict_alloc(dict->n, dict->n, CM_DG_DICT_IDENTITY, NULL);
    cm_dg_fill_identity(id_dict->data, id_dict->n);

    printf("Identity dictionary.\n");
    /* Print dictionary out */
    for (i = 0; i < id_dict->n; i++) {
        for (j = 0; j < id_dict->p; j++)
            printf("%1.3f ", cm_dg_get_c(id_dict, i, j)); 
        printf("\n");
    }

    /* Make a combined dictionary */
    cm_dg_dict_t *comb_dict = cm_dg_dict_combine(dict, id_dict, NULL);

    /* Free old dictionaries */
    cm_dg_dict_free(dict);
    cm_dg_dict_free(id_dict);

    printf("Combined dictionary.\n");
    /* Print dictionary out */
    for (i = 0; i < comb_dict->n; i++) {
        for (j = 0; j < comb_dict->p; j++)
            printf("%1.3f ", cm_dg_get_c(comb_dict, i, j)); 
        printf("\n");
    }

    cm_dg_dict_free(comb_dict);
    return(0);

}

