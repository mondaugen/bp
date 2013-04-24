/* Copyright 2013 Nicholas Esterer. All Rights Reserved. */
/* Routines for generating waveform dictionaries for algorithms like Basis
 * Pursuit, Matching Pursuit, etc. */

#include "cm_dict_gen.h" 
#include <math.h>
#include <assert.h> 
#include <stdlib.h>
#include <dirent.h>
#include <sys/types.h> 
#include <string.h> 
#include <sndfile.h> 
#include <fnmatch.h> 

#ifdef CM_DEBUG
 #include <stdio.h> 
#endif

/* transpose a matrix in-place, taken from:
 * http://rosettacode.org/wiki/Matrix_transposition#C
 */
void
cm_dg_matrix_transpose(double *m, int w, int h)
{
	int start, next, i;
	double tmp;
 
	for (start = 0; start <= w * h - 1; start++) {
		next = start;
		i = 0;
		do {	i++;
			next = (next % h) * w + next / h;
		} while (next > start);
		if (next < start || i == 1) continue;
 
		tmp = m[next = start];
		do {
			i = (next % h) * w + next / h;
			m[next] = (i == start) ? tmp : m[i];
			next = i;
		} while (next > start);
	}
}

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
    memset(result->data, 0, sizeof(double)*n*p);
    result->n = n;
    result->p = p;
    result->type = type;
    result->info = info;
    return result;
}

/* Combine two dictionaries to form a new dictionary */
/* a and b are two dictionaries and infocombine is a function that will combine
 * the two info fields or NULL if you don't care about the info field in the
 * two dictionaries */
cm_dg_dict_t *
cm_dg_dict_combine(cm_dg_dict_t *a, cm_dg_dict_t *b,
        void* (*infocombine)(void*,void*))
{
    /* The lengths of the signals in the two dictionaries must be equal */
    if (a->n != b->n)
        return NULL;
    cm_dg_dict_t *result;
    assert((result = cm_dg_dict_alloc(a->n, a->p + b->p, CM_DG_DICT_CUSTOM,
                    (infocombine == NULL) ? NULL : infocombine(a->info,
                        b->info))));
    size_t i;
    for (i = 0; i < result->n; i++) {
        memcpy( (result->data + i*(a->p + b->p)),
                (a->data + i*(a->p)),
                a->p * sizeof(double) );
        memcpy( (result->data + i*(a->p + b->p) + a->p),
                (b->data + i*(b->p)),
                b->p * sizeof(double));
    }
    return result;
}

/* Count the number of files in the directory dir that match the wildcard
 * pattern wc */
size_t
cm_dg_wc_count(char *dir, char *wc)
{
    DIR *dp;
    struct dirent *ep;
    size_t count = 0;

    if ((dp = opendir(dir)) != NULL)
        while (ep = readdir(dp))
            if (fnmatch(wc, ep->d_name, FNM_CASEFOLD) != FNM_NOMATCH)
                count++;

    return count;
}

/* Allocate a dictionary of of waveforms loaded from the directory dir that
 * match the wildcard pattern wc.
 * n is the number of samples to read from the waveforms. If a file contains
 * less than this number of samples, then the remaining samples will be 0.
 * info is set to NULL */
cm_dg_dict_t *
cm_dg_dict_alloc_from_dir(size_t n, char *dir, char *wc)
{
    DIR *dp;
    struct dirent *ep;

    size_t count = cm_dg_wc_count(dir, wc);

#ifdef CM_DEBUG
    fprintf(stderr,"Count = %lu\n", count);
#endif  

    cm_dg_dict_t *result;
    assert((result = cm_dg_dict_alloc(n, count, CM_DG_DICT_CUSTOM, NULL))
            != NULL);

    SNDFILE *sndfile;
    SF_INFO sfinfo;

    if ((dp = opendir(dir)) != NULL) {
        size_t i = 0;
        while (ep = readdir(dp)) {
            if ((fnmatch(wc, ep->d_name, FNM_CASEFOLD) != FNM_NOMATCH)
                   && (i < count)) {
                i++;
                char fullpath[strlen(dir) + strlen(ep->d_name) + 2];
                strcpy(fullpath,dir);
                strcat(fullpath,"/");
                strcat(fullpath,ep->d_name);
                memset(&sfinfo, 0, sizeof(SF_INFO));
                if (!(sndfile = sf_open(fullpath, SFM_READ, &sfinfo))) {
                    fprintf(stderr,"Couldn't open %s, skipping.\n",
                            fullpath);
                    continue;
                }
                if (sfinfo.channels != 1){
                    fprintf(stderr,"Number of channels of %s is %d and not "
                            "1, skipping.\n", fullpath, sfinfo.channels);
                    continue;
                }
#ifdef CM_DEBUG
                fprintf(stderr, "Reading file %s\n", fullpath);
#endif  
                /* Read the soundfiles into the rows */
                sf_read_double(sndfile, (result->data) + (i-1)*n ,n);
                sf_close(sndfile);
            }
        }
    }
    /* swap rows and columns so the dictionary is addressed correctly */
    cm_dg_matrix_transpose(result->data, n, count);
#ifdef CM_DEBUG 
    size_t i, j;
    for (i = 0; i < result->n; i++) {
        for (j = 0; j < result->p; j++)
            fprintf(stderr,"%1.3f ", cm_dg_get_c(result, i, j)); 
        fprintf(stderr,"\n");
    }
#endif  
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

