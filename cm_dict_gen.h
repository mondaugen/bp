#ifndef CM_DICT_GEN_H
#define CM_DICT_GEN_H 

#include <stddef.h> 

typedef enum __cm_dg_dict_type {
    CM_DG_DICT_FOURIER,
    CM_DG_DICT_IDENTITY,
    CM_DG_DICT_CUSTOM
} cm_dg_dict_type;

typedef struct __cm_dg_dict {
    size_t n;   /* length of each waveform in dictionary */
    size_t p;   /* number of waveforms in dictionary */
    cm_dg_dict_type type;   /* type of dictionary */
    void *info; /* pointer to some extra info about the dictionary */
    double data[0]; /* waveforms of dictionary are stored here */
} cm_dg_dict_t;

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
cm_dg_fill_fourier(double *phi, size_t n, size_t l);

/* Fills the dictionary with elements of an n x n identity matrix. */
int
cm_dg_fill_identity(double *phi, size_t n);

/* Allocate a dictionary and return a pointer to it. info must be a pointer to
 * freeable memory (memory allocated on the heap), or NULL */
cm_dg_dict_t *
cm_dg_dict_alloc(size_t n, size_t p, int type, void *info);

void
cm_dg_dict_free(cm_dg_dict_t *d);

/* get the data in the i-th row and j-th column with c-style indexing */
double
cm_dg_get_c(cm_dg_dict_t *d, size_t i, size_t j);

/* get the data in the i-th row and j-th column with mathematics indexing
 * (starting from 1) */
double
cm_dg_get_m(cm_dg_dict_t *d, size_t i, size_t j);

/* Allocate a dictionary of of waveforms loaded from the directory dir that
 * match the wildcard pattern wc.
 * n is the number of samples to read from the waveforms. If a file contains
 * less than this number of samples, then the remaining samples will be 0.
 * info is set to NULL */
cm_dg_dict_t *
cm_dg_dict_alloc_from_dir(size_t n, char *dir, char *wc);

/* Combine two dictionaries to form a new dictionary */
/* a and b are two dictionaries and infocombine is a function that will combine
 * the two info fields or NULL if you don't care about the info field in the
 * two dictionaries */
cm_dg_dict_t *
cm_dg_dict_combine(cm_dg_dict_t *a, cm_dg_dict_t *b,
        void* (*infocombine)(void*,void*));

#endif /* CM_DICT_GEN_H */
