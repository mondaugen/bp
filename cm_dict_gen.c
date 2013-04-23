/* Copyright 2013 Nicholas Esterer. All Rights Reserved. */
/* Routines for generating waveform dictionaries for algorithms like Basis
 * Pursuit, Matching Pursuit, etc. */

#include "cm_dict_gen.h" 
#include <math.h>
#include <assert.h> 
#include <stdlib.h> 

#ifdef CM_DEBUG
 #include <stdio.h> 
#endif


/* Fill an array of doubles with an overcomplete Fourier dictionary.
 * phi is an array of doubles in which to store the dictionary.
 * phi must have size >= (l*n)*n*sizeof(double)
 * phi represents a matrix and can be used with the cm_dg_dict struct to ease
 * access of elements.
 * The elements of phi are stored as follows:
 * The l*n waveforms make up the columns (of length n) so the ith sample of the
 * jth waveform is found at phi[i*l*n + j].
 * Returns non-zero on error.
 */
int
cm_dg_fill_fourier(double *phi, size_t n, size_t l)
{
    size_t i, j; /* rows and columns, resp. */
    double omega;
    
    for (j = 0; j <= ((l*n)/2); j++) {
    
        /* fill the first (l*n)/2 + 1 columns with cosines */
        omega = 2. * M_PI * (double)j / ((double)n*(double)l);

#ifdef CM_DEBUG
       fprintf(stderr,"Omega = %f at j = %lu\n", omega, j);
#endif  
    
        for (i = 0; i < n; i++) {
            phi[i*l*n + j] = cos(omega * (double)i);
#ifdef CM_DEBUG
           fprintf(stderr,"phi(%lu,%lu) = %f\n", i, j, phi[i*l*n + j]);
#endif  
        }
    
    }

    for ( ; j < (l*n); j++ ) {

        /* fill column (l*n)/2 + 1 to (l*n) with sines */
        omega = 2. * M_PI * (double)(j - ((l*n)/2)) / ((double)n*(double)l);
#ifdef CM_DEBUG
       fprintf(stderr,"Omega = %f at j = %lu\n", omega, j);
#endif  

        for (i = 0; i < n; i++) {
            phi[i*l*n + j] = sin(omega * (double)i);
#ifdef CM_DEBUG
           fprintf(stderr,"phi(%lu,%lu) = %f\n", i, j, phi[i*l*n + j]);
#endif  
        }
        
    }

    return(0);

}

/* Fills the dictionary with elements of an n x n identity matrix. */
int
cm_dg_fill_identity(double *phi, size_t n)
{
    size_t i, j;
    for (j = 0; j < n; j++) {
	for (i = 0; i < n; i++) {
	    if (i == j)
		phi[i*n + j] = 1;
	    else
		phi[i*n + j] = 0;
	}
    }
    return(0);
}


/* Allocate a dictionary and return a pointer to it. info must be a pointer to
 * freeable memory (memory allocated on the heap), or NULL */
cm_dg_dict_t *
cm_dg_dict_alloc(size_t n, size_t p, int type, void *info)
{
    cm_dg_dict_t *result;
    assert((result = (cm_dg_dict_t*)malloc(sizeof(cm_dg_dict_t)
            + sizeof(double)*n*p)));
    result->n = n;
    result->p = p;
    result->type = type;
    result->info = info;
    return result;
}

void
cm_dg_dict_free(cm_dg_dict_t *d)
{
    if(d->info != NULL)
        free(d->info);
    free(d);
}

/* get the data in the i-th row and j-th column with c-style indexing */
double
cm_dg_get_c(cm_dg_dict_t *d, size_t i, size_t j)
{
    /* check bounds */
    assert((i < d->n) && (j < d->p));
    return d->data[i*(d->p) + j];
}

/* get the data in the i-th row and j-th column with mathematics indexing
 * (starting from 1) */
double
cm_dg_get_m(cm_dg_dict_t *d, size_t i, size_t j)
{
    /*check bounds*/
    assert((i <= d->n) && (j <= d->p));
    return d->data[(i-1)*(d->p) + j - 1];
}

