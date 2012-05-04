#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "matrix.h"
#include "linalg.h"


matrix *matrix_inverse(matrix *A) {

	matrix *R;
	matrix *invA;
	if (A->n != A->m) {
		fprintf(stderr, "Cannot take inverse of non-square matrix\n");
		exit(-1);
	}
	R = QR_decomposition(A);
	invA = matrix_inverse_QR(A,R);
	matrix_free(R);
	// If we implement an in-place matrix transposition
	// algorithm, the following should be avoided.
	// Actually, here we simulate an in-place algorithm
	free(A->root);
	free(A->a);
	A->root = invA->root;
	A->a = invA->a;
	// n and m remain unchanged
	free(invA);
	return A;
}

matrix *matrix_inverse_QR(matrix *Q, matrix *R) {

	if (Q->n != Q->m) {
		fprintf(stderr, "Cannot take inverse of non-square matrix.\n");
		exit(-1);
	}

	return matrix_inverse_QR_prealloc(matrix_create(Q->n,Q->n), Q, R);
}

matrix *matrix_inverse_prealloc(matrix *invA, matrix *A) {

	matrix *R;

	if (A->n != A->m) {
		fprintf(stderr, "Cannot take inverse of non-square matrix\n");
		exit(-1);
	} else if (invA->n != A->n || invA->m != A->n) {
		fprintf(stderr, "Wrong dimensions on preallocated inverse matrix.\n");
		exit(-1);
	}

	R = QR_decomposition(A);
	matrix_inverse_QR_prealloc(invA, A, R);
	free(R);

	return invA;
}

matrix *matrix_inverse_QR_prealloc(matrix *invA, matrix *Q, matrix *R) {

	size_t i, j, k;

	if (Q->n != Q->m) {
		fprintf(stderr, "Cannot take inverse of non-square matrix.\n");
		exit(-1);
	} else if (invA->n != Q->n || invA->m != Q->n) {
		fprintf(stderr, "Wrong dimensions on preallocated inverse matrix.\n");
		exit(-1);
	}


	// The following is a cache-optimized version of the linalg-solve
	// routine that solves for all columns in a matrix instead
	// of just 1 column. This way, we can do each row 
	// before the next iteration in back-substitution
	// to enhance locality

	// The input columns form the identity matrix,
	// so we just need to solve Rx = Q^T[:,col] colum by column

	//	x[:][i] = Q^T[:][i] = Q[i][:]

	// If we had an in-place matrix transposition
	// algorithm (like the cache-oblivious one described on Wikipedia),
	// we could probably make this faster
	// and use up 1/3 less memory.

	i = R->m;

	do {  // Loop over rows

		i--;	
		for (k = 0; k < Q->n; k++)
			invA->a[i][k] = Q->a[k][i];
	
		for (j = R->m - 1; j > i; j--) // Backsub loop
			for (k = 0; k < Q->n; k++) // Column
			   invA->a[i][k] -= R->a[i][j] * invA->a[j][k];

		for (k = 0; k < Q->n; k++) // Column
			invA->a[i][k] /= R->a[i][i];
	} while (i != 0);
	
	return invA;
}


