#ifndef CM_BP_H
#define CM_BP_H 

#include <stddef.h> 

#define CM_BP_LP_VERBOSE 0x1

/* Given a dictionary D, of n rows and p columns (n is signal length, p is
 * number of items in the dictionary), a signal s of length n and a coefficient
 * vector x of length p, solve the linear program min|x| subject to Dx = s.
 * D is stored in row-major order. */
int
cm_bp_lp_solve(double *D, size_t n, size_t p, double *x, double *s,
        unsigned int solveflags);

#endif /* CM_BP_H */
