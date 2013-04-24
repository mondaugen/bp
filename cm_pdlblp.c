/* An implementation of the "primal-dual log-barrier linear programming
 * algorithm" of Chen, Donoho and Saunders (1998)
 *
 * Nicholas Esterer, April 2013, for MUMT 622.
 */

#include <cm_pdlblp.h>
#include <math.h>
#include <stdio.h> 

/* Multiply each row of matrix M by an element of vector v
 * If V is a matrix with v on the main diagonal, then the result is M = VM */
int
_matrix_rows_mul(gsl_matrix *M, gsl_vector *v)
{
    if (M->size1 != v->size)
        return(-1); /* number of rows must equal vector size */
    gsl_vector_view row_view;
    size_t i;
    for (i = 0; i < M->size1; i++) {
        row_view = gsl_matrix_row(M,i);
        gsl_vector_scale(&row_view.vector, gsl_vector_get(v,i));
    }
    return(0);
}

/* find max{ p : x + p*dx >= 0 } */
int
_p_max_find(double *p, gsl_vector *x, gsl_vector *dx)
{
    double lob = DBL_MIN, hib = DBL_MAX; /* lower and upper bounds on p*/
    size_t i;
    
    for (i = 0; i < x->size; i++) {
        if (gsl_vector_get(dx,i) == 0){ /* There is no bound on p if x >= 0 */
            if (gsl_vector_get(x,i) < 0)
                return(-1); /* No such p exists */
            continue;
        }
        double ptmp = -1. * gsl_vector_get(x,i) / gsl_vector_get(dx,i);
        if (gsl_vector_get(dx,i) < 0) /* There is an upper bound on p */
            if (ptmp < hib)
                hib = ptmp;
        if (gsl_vector_get(dx,i) > 0) /* There is an lower bound on p */
            if (ptmp > lob)
                lob = ptmp;
        if (hib < lob) /* If upper bound lower than lower bound, no p exists */
            return(-2);
    }
    *p = hib; /* Result goes in p */
    return(0);
}

/* calculate t. t, c, x and z have length p (dictionary length), y is length n
 * (signal length) and A is size (n,p) */
int
_t_calc(gsl_vector *t, gsl_vector *c, double gamma, gsl_vector *x, gsl_vector *z,
        gsl_matrix *A, gsl_vector *y)
{
    /* t = c */
    gsl_vector_memcpy(t, c);
    /* t = gamma^2*x + c */
    gsl_blas_daxpy(gamma*gamma, x, t);
    /* t = gamma^2*x + c - z */
    gsl_blas_daxpy(-1, z, t);
    /* t = -A^T*y + t gamma^2*x + c - z */
    gsl_blas_dgemv(CblasTrans, -1, A, y, 1., t);
    return(0);
}

/* calculate r. r, b and y size n, x size p, A size (n,p) */
_r_calc(gsl_vector *r, gsl_vector *b, gsl_matrix *A, gsl_vector *x,
        double delta, gsl_vector *y)
{
    /* r = b */
    gsl_vector_memcpy(r, b);
    /* r = b - delta^2*y */
    gsl_blas_daxpy(-1.*delta*delta, y, r);
    /* r = b - delta^2*y - Ax */
    gsl_blas_dgemv(CblasNoTrans, -1, A, x, 1., r);
    return(0);
}

/* calculate v. v, z and x size p */
_v_calc(gsl_vector *v, double mu, gsl_vector *z, gsl_vector *x)
{
    /* v = z */
    gsl_vector_memcpy(v, z);
    /* v = Z*x */
    gsl_vector_mul(v,x);
    /* v = -Zx */
    gsl_vector_scale(v,-1.);
    /* v = mu*e - Zx */
    gsl_vector_add_constant(v, mu);
    return(0);
}

/* D is a diagonal matrix size (p,p). we store D as a vector, length p, to save
 * space, x and z are vectors of size p */
_D_calc(gsl_vector *D, gsl_vector *x, gsl_vector *z, double gamma)
{
    size_t i;
    for (i = 0; i < D->size; i++)
        gsl_vector_set(D, i, 1. / (gsl_vector_get(z,i) / gsl_vector_get(x,i)
                    + gamma*gamma));
    return 0;
}

/* calculate delta_y (dy). dy is a vector of size n, D is a vector of size p, A
 * is size (n,p), t is size p, x is size p, v is size p, r is size n, delta is a
 * double, B is a matrix size (p+n,n), V is matrix size (n,n), g is vector size
 * (n+p) s is vector size(n), work is vector size (n) */
_delta_y_calc(gsl_vector *dy, gsl_vector *D, gsl_matrix *A, gsl_vector *t,
        gsl_vector *x, gsl_vector *v, gsl_vector *r, double delta, gsl_matrix *B,
        gsl_matrix *V, gsl_vector *g, gsl_vector *s, gsl_vector *work)
{
    size_t i, n, p;
    n = A->size1; /* The number of rows in A */
    p = A->size2; /* The number of columns in A */
    
    /* Split matrix B into top and bottom views */
    gsl_matrix_view B_top = gsl_matrix_submatrix(B,0,0,p,n);
    gsl_matrix_view B_bottom = gsl_matrix_submatrix(B,p,0,n,n);
    
    /* Copy transposed A into B */
    gsl_matrix_transpose_memcpy(&B_top.matrix, A);

    /* Put D^(1/2) into the top of g */
    for (i = 0; i < p; i++)
        gsl_vector_set(g,i,sqrt(gsl_vector_get(D,i)));

    /* Get view of top of g only */
    gsl_vector_view g_view = gsl_vector_subvector(g, 0, p);
    
    /* set B_top = D^(1/2)B_top, or B_top = D^(1/2)A^T */
    _matrix_rows_mul(&B_top.matrix, &g_view.vector);
    
    /* set B_bottom = delta*I */
    gsl_matrix_set_identity(&B_bottom.matrix);
    gsl_matrix_scale(&B_bottom.matrix, delta);

    /* set g_top = g_top * (t - X^(-1)*v) */
    for (i = 0; i < p; i++)
        gsl_vector_set(g, i, gsl_vector_get(g,i)*(gsl_vector_get(t,i)
                    - gsl_vector_get(v,i) / gsl_vector_get(x,i)));

    /* set g_bottom = r / delta */
    for ( ; i < p+n; i++)
        gsl_vector_set(g, i, gsl_vector_get(r,i - p) / delta);

    /* solve B*dy = g in the least squares sense */
    gsl_linalg_SV_decomp(B,V,s,work);
    gsl_linalg_SV_solve(B,V,s,g,dy);

    /* answer is now in dy */
    return(0);
}

/* calculate dx */
_dx_calc(gsl_vector *dx, gsl_vector *D, gsl_matrix *A,
        gsl_vector *dy, gsl_vector *x, gsl_vector *v, gsl_vector *t)
{
    /* dx = X^(-1)v - t */
    gsl_vector_memcpy(dx,v);
    gsl_vector_div(dx,x);
    gsl_vector_sub(dx,t);
    /* dx = A^T*dy + X^(-1)v - t */
    gsl_blas_dgemv(CblasTrans, 1, A, dy, 1, dx);
    /* dx = D*A^T*dy + D*(X^(-1)v - t) */
    gsl_vector_mul(dx,D);
    return(0);
}

/* calculate dz */
_dz_calc(gsl_vector *dz, gsl_vector *x, gsl_vector *v, gsl_vector *z,
        gsl_vector *dx)
{
    /* dz = dx */
    gsl_vector_memcpy(dz, dx);
    /* dz = Z * dx */
    gsl_vector_mul(dz, z);
    /* dz = -Z * dx */
    gsl_vector_scale(dz,-1.);
    /* dz = v - Z * dx */
    gsl_vector_add(dz, v);
    /* dz = X^(-1)v - X^(-1)*Z*dx */
    gsl_vector_div(dz, x);
    return(0);
}


/* solve min{c^Tx + 1/2*|gamma*x|^2 + 1/2|p|^2} subject to Ax + delta*p = b,
 * x >= 0.
 * For BP, c is usually e (a vector of all ones), A(n,p) is the dictionary of p
 * waveforms of length n, b(n) is the signal and the resulting x(p) is the vector of
 * coefficients characterizing the basis. Gamma and delta are the regularization
 * parameters "feasibility tolerance" and "duality gap tolerance", K is the
 * maximum number of iterations of the algorithm to carry out.
 */
int
pdlblp_solve(gsl_vector *x, gsl_matrix *A, gsl_vector *b, gsl_vector *c,
        double gamma, double delta, size_t K, double yinit,
	double zinit, double muinit)
{
    size_t p, n, k;
    double pp, pd, ztx, mu, pig, dig, dg;
    int err = 0;

    n = A->size1;
    p = A->size2;

    gsl_matrix *B, *V;
    B = gsl_matrix_alloc(p+n,n);
    V = gsl_matrix_alloc(n,n);

    gsl_vector *t, *y, *z, *r, *v, *D, *dx, *dy, *dz, *g, *s, *work;

    t    = gsl_vector_alloc(p);
    y    = gsl_vector_alloc(n);
    z    = gsl_vector_alloc(p);
    r    = gsl_vector_alloc(n);
    v    = gsl_vector_alloc(p);
    D    = gsl_vector_alloc(p);
    dx   = gsl_vector_alloc(p);
    dy   = gsl_vector_alloc(n);
    dz   = gsl_vector_alloc(p);
    g    = gsl_vector_alloc(p+n);
    s    = gsl_vector_alloc(n);
    work = gsl_vector_alloc(n);

    /* initialize x, y, z, mu */
//    gsl_vector_set_all(x,xinit);
    gsl_vector_set_all(y,yinit);
    gsl_vector_set_all(z,zinit);
    mu = muinit;

    for(k = 0; k < K; k++){

	fprintf(stderr,"In [%s]: k = %lu\n",__FILE__,k);

        /* calculate t */
	_t_calc(t, c, gamma, x, z, A, y);

        /* calculate r */
        _r_calc(r, b, A, x, delta, y);

        /* calculate v */
        _v_calc(r, mu, z, x);

        /* calculate D */
        _D_calc(D, x, z, gamma);

        /* calculate delta y */
        _delta_y_calc(dy, D, A, t, x, v, r, delta, B, V, g, s, work);

        /* calculate dx */
        _dx_calc(dx, D, A, dy, x, v, t);

        /* calculate dz */
        _dz_calc(dz, x, v, z, dx);

        /* calculate primal step size pp */
        if((err = _p_max_find(&pp, x, dx)))
	    fprintf(stderr,"Error %d at [%s]:%d\n",err,__FILE__,__LINE__);
        pp *= .99;

        /* calculate dual step size pd */
        if((err = _p_max_find(&pd, z, dz)))
	    fprintf(stderr,"Error %d at [%s]:%d\n",err,__FILE__,__LINE__);
        pd *= .99;

        /* update x */
        gsl_blas_daxpy(pp,dx,x);

        /* update y */
        gsl_blas_daxpy(pd,dy,y);

        /* udpate z */
        gsl_blas_daxpy(pd,dz,z);

        /* update mu = (1 - min(pp, pd, .99))*mu */
        mu *= (1. - GSL_MIN( GSL_MIN( pp, pd ), .99 ));

        /* calc z^T*x */
        gsl_blas_ddot(x, z, &ztx);

	pig = (gsl_blas_dnrm2(r) / (1. + gsl_blas_dnrm2(x)));
	fprintf(stderr,"Primal infeasibility gap: %f.\n", pig);

	dig = (gsl_blas_dnrm2(t) / (1. + gsl_blas_dnrm2(y)));
	fprintf(stderr,"Dual infeasibility gap: %f.\n", dig);

	dg  = (ztx / (1. + gsl_blas_dnrm2(z)*gsl_blas_dnrm2(x)));
	fprintf(stderr,"Duality gap: %f.\n", dg);

        /* check primal/dual infeasibility and duality gap */
        if (( pig < gamma)
                && ( dig < gamma )
                && ( dg < delta ))
            goto cleanup; /* suitable x has been found */

    }

cleanup:

    gsl_matrix_free(B);
    gsl_matrix_free(V);
    gsl_vector_free(t);
    gsl_vector_free(y);
    gsl_vector_free(z);
    gsl_vector_free(r);
    gsl_vector_free(v);
    gsl_vector_free(D);
    gsl_vector_free(dx);
    gsl_vector_free(dy);
    gsl_vector_free(dz);
    gsl_vector_free(g);
    gsl_vector_free(s);
    gsl_vector_free(work);

    return(err);
}
