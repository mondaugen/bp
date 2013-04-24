/* Test the loading of soundfiles into a dictionary */

#include <cm/sndfio.h> 
#include <cm_dict_gen.h>
#include <stdlib.h> 
#include <stdio.h>

int
main(void)
{
    /* write some junk to some sound-files */
    size_t i,j;
    double junk[10];
    char filename[] = "./a.wav";
    for (i = 0; i < 5; i++) {
        /* 5 soundfiles */
        for ( j = 0; j < 10; j++)
            junk[j] = ((double)(i*10 + j)) * 0.01;
        cm_write_sf_from_dbl_ary_params(junk, filename, 5, 44100, 1,
                (SF_FORMAT_WAV | SF_FORMAT_PCM_16));
        filename[2] += 1;
    }

    /* Load the junk from the soundfiles */
    cm_dg_dict_t *dict = cm_dg_dict_alloc_from_dir(10,"./","*.wav");

    for (i = 0; i < dict->n; i++) {
        for (j = 0; j < dict->p; j++)
            printf("%1.3f ", cm_dg_get_c(dict, i, j)); 
        printf("\n");
    }

    cm_dg_dict_free(dict);
    return(0);

}

