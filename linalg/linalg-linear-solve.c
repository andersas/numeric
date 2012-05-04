#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"
#include "linalg.h"


// Solve the linear equation Ax = b where QR = A is the QR decomposition
// of A and b and x are column vectors. Caller should supply x to be filled.

void solve_linear_system(matrix *Q, matrix *R, matrix *b, matrix *x) {

	double sum;
	size_t i,j;

	if (b->m != 1 || x->m != 1 || x->n != R->m || Q->n != b->n) {
		fprintf(stderr, "Wrong dimensions on input vectors to linear solver.\n");	
		exit(-1);
	}

	// Solving Ax=b is equivalent to solving 
	// Rx = Q^T*b.
	
	// Store Q^T * b in x[]

	for (i = 0; i < R->m; i++) {
		x->a[i][0] = matrix_dot_cols(Q, b, i, 0);
	}
	// Solve the system by back substitution
	
	i = R->m;
	do {
		i--;
		sum = 0.0;
		for (j = R->m - 1; j > i; j--) {
			sum += R->a[i][j] * x->a[j][0];
		}
		x->a[i][0] = (x->a[i][0] - sum) / R->a[i][i];

	} while (i != 0);

}


