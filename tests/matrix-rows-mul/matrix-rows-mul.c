/* Testing matrix rows multiply routine */

#include <cm_pdlblp.h> 
#include <stdio.h>

int
main(void)
{
    double A[] = {  1, 2, 3,
                    4, 5, 6,
                    7, 8, 9 };
    gsl_matrix_view A_view = gsl_matrix_view_array(A, 3, 3);
    double v[] = {  3, 2, 1 };
    gsl_vector_view v_view = gsl_vector_view_array(v, 3);

    _matrix_rows_mul(&A_view.matrix, &v_view.vector);
    gsl_matrix_fprintf(stdout, &A_view.matrix, "%f");

    return(0);
}
