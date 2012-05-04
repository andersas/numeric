#include <stdio.h>
#include <stdlib.h>
#include "../linalg/linalg.h"
#include "leastSquares.h"


// Makes the linear least squares solution to the system
// (5) in the notes.
// covMat is a preallocated mxm matrix, m i
//

void leastSquares(size_t n, size_t m, matrix *x, matrix *y, matrix *dy, double (*f[m])(double), matrix *coeff, matrix *cov) {

	matrix *A, *b;
	size_t i,k;

	if (x->n != n || x->m != 1 || x->n != y->n || x->m != y->m || dy->n != x->n || dy->m != x->m) {
		fprintf(stderr,"leastSquares: Wrong matrix dimensions!\n");
		exit(-1);
	}

	if (coeff->n != m || coeff->m != 1 || cov->n != m || cov->m != m) {
		fprintf(stderr,"leastSquares: Wrong output matrix dimensions!\n");
		exit(-1);
	}

	if (n < m) {
		fprintf(stderr,"leastSquares: Insufficient data for fit.\n");
	}

	A = matrix_create(n,m);
	b = matrix_create(n,1);

	// Fill out a matrix of function values for solving
	
	for (i = 0; i < n; i++) {
		for (k = 0; k < m; k++) {
			A->a[i][k] = f[k](x->a[i][0])/dy->a[i][0];
		}
		b->a[i][0] = y->a[i][0]/dy->a[i][0];
	}

	// Then do a QR decomposition of A

	QR_decomposition_prealloc(A,cov);

	// and backsubstitute
	
	solve_linear_system(A,cov,b,coeff);

	// The covariance matrix is calculated as (R^T * R)^-1.
	// R is currently stored in cov.
	
	matrix_free(A);	
	matrix_free(b);	

	A = matrix_create(m,m);

	for (i = 0; i < m; i++)
		for (k = 0; k < m; k++)
			A->a[i][k] = matrix_dot_cols(cov,cov,i,k);

	matrix_inverse_prealloc(cov, A);
	matrix_free(A);

}

