#ifndef LINALG
#define LINALG

#include "matrix.h"
#include "pp.h"

matrix * QR_decomposition(matrix *A);
matrix * QR_decomposition_prealloc(matrix *A, matrix *R);

void solve_linear_system(matrix *Q, matrix *R, matrix * b, matrix * x);

double QR_absdet(matrix *Q, matrix *R);

matrix *matrix_inverse(matrix *A);
matrix *matrix_inverse_QR(matrix *Q, matrix *R);
matrix *matrix_inverse_prealloc(matrix *invA, matrix *A);
matrix *matrix_inverse_QR_prealloc(matrix *invA, matrix *Q, matrix *R);



#endif
