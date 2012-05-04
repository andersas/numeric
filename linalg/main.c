#include <stdio.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "matrix.h"
#include "linalg.h"
#include "pp.h"

int main(int argc, char **argv) {

	matrix *A, *B, *Q, *R, *b, *x;
	gsl_matrix_view m;
	gsl_vector *tau;


if (argc != 1) {
	Q = matrix_create_random(3000,3000);
	
	if (!strcmp(argv[1],"gsl")) {
		m = gsl_matrix_view_array(Q->root, Q->n, Q->m);
		tau = gsl_vector_alloc (Q->n);
		gsl_linalg_QR_decomp (&m.matrix, tau);
		
	} else {
		R = QR_decomposition(Q);
	}

} else {
printf("Detected %zu logical processors.\n", numCPUs());

	printf("Loading A from A.inp...\n");
	A = matrix_load("A.inp");
	Q = matrix_copy(A);
	printf("Loading b from b.inp...\n");
	b = matrix_load("b.inp");

	printf("Doing QR decomposition on A...\n");
	R = QR_decomposition(Q);

	printf("Calculating Q x R...\n");
	B = matrix_multiply(Q,R);
	printf("Writing Q.out, R.out, A.out (=Q x R)...\n");
	matrix_store(Q, "Q.out");
	matrix_store(R, "R.out");
	matrix_store(B, "A.out");

	matrix_free(B);

	printf("Solving Ax = b...\n");
	x = matrix_create(A->n,1);
	printf("doing it..\n");
	solve_linear_system(Q,R,b,x);
	printf("Outputting result to x.out...\n");
	matrix_store(x, "x.out");
	matrix_free(x);
	matrix_free(b);

	printf("abs(det(A)) = %f\n", QR_absdet(Q,R));

	printf("Calculating the inverse of A...\n");
	B = matrix_inverse_QR(Q,R);
	printf("Outputting inverse of A to invA.out...\n");
	matrix_store(B,"invA.out");

	matrix_free(B);
	matrix_free(Q);
	matrix_free(R);
	matrix_free(A);
}

	return 0;

}
