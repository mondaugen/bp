/* Testing _p_max_find routine */

#include <gsl/gsl_vector.h>
#include <cm_pdlblp.h>
#include <stdio.h> 

int
main(void)
{
    gsl_vector *x = gsl_vector_alloc(4);
    gsl_vector *dx = gsl_vector_alloc(4);

    gsl_vector_set(x,0,-4);
    gsl_vector_set(x,1,-5);
    gsl_vector_set(x,2,-1);
    gsl_vector_set(x,3,-1);

    gsl_vector_set(dx,0,2);
    gsl_vector_set(dx,1,-2);
    gsl_vector_set(dx,2,0.5);
    gsl_vector_set(dx,3,-0.5);

    double p;
    int err;
    if ((err = _p_max_find(&p, x, dx)))
        fprintf(stderr,"Error %d in finding p.\n", err);
    else
        fprintf(stderr,"p is %f \n", p);

    gsl_vector_scale(x,-1);
    gsl_vector_scale(dx,-1);

    if ((err = _p_max_find(&p, x, dx)))
        fprintf(stderr,"Error %d in finding p.\n", err);
    else
        fprintf(stderr,"p is %f \n", p);

    gsl_vector_free(x);
    gsl_vector_free(dx);

    return(0);
}
