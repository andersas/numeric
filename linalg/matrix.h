#ifndef MATRIX
#define MATRIX
#include <stdlib.h>

typedef struct matrix {

	size_t n,m;
	double *root; // In C++, this would be private.
	double **a;

} matrix;

double matrix_dot_product(matrix *A, size_t row, matrix *B, size_t col);
double matrix_single_dot_product(size_t len, double *v, double *w);
double matrix_dot_cols(matrix *A, matrix *B, size_t colA, size_t colB);
double matrix_vector_norm(matrix *v);
matrix *matrix_create(size_t n, size_t m);
matrix *matrix_create_zero(size_t n, size_t m);
matrix *matrix_create_random(size_t n, size_t m);
matrix *matrix_add(matrix *A, matrix *B);
matrix *matrix_sub(matrix *A, matrix *B);
matrix *matrix_scale(matrix *A, double b);
matrix *matrix_scale_col(matrix *A, size_t col, double b);
matrix *matrix_scale_row(matrix *A, size_t row, double b);
matrix *matrix_create_diagonal(size_t n, matrix *diagonal);
matrix *matrix_create_diagonal_scalar(size_t n, double diagonal);
matrix *matrix_identity(size_t n);
matrix *matrix_get_diagonal(matrix *mat, matrix *storage);
matrix *matrix_set_diagonal(matrix *mat, const matrix *storage);
matrix *matrix_set_diagonal_scalar(matrix *mat, double diagonal);
matrix *matrix_zero(matrix *mat);
matrix *matrix_swap_rows(matrix *A, size_t i, size_t j);
matrix *matrix_swap_cols(matrix *A, size_t i, size_t j);
void matrix_free(matrix * mat);
matrix *matrix_column_from_array(size_t len, const double v[]);
matrix *matrix_row_from_array(size_t len, const double v[]);
matrix *matrix_clone(matrix *mat);
matrix *matrix_copy(matrix *mat);
matrix *matrix_clone_preallocated(matrix *dest, matrix *src);
matrix *matrix_copy_preallocated(matrix *dest, matrix *src);
matrix *matrix_store(matrix *mat, char *filename);
matrix * matrix_load(char *filename);

// More advanced functions:

matrix *matrix_transpose(matrix *mat);

matrix *matrix_multiply(matrix *A, matrix *B);
matrix *matrix_multiply_prealloc(matrix *C, matrix *A, matrix *B);


#endif
