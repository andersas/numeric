#include <stdio.h>
#include "../linalg/linalg.h"


#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

// Helper function to create a square jacobian matrix
matrix *square_jacobian(size_t n, double (*f[n])(matrix *x), matrix *J, matrix *x, matrix *minusfx, double dx) {
size_t i, j;
// Calculate the jacobian of the square system f(x) in x.
// J is an output matrix containing the jacobian


// The jacobian is a column vector of gradients of the coordinate functions.	
	for (i = 0; i < n; i++) // loop over rows
		for (j = 0; j < n; j++) { // Loop over columns
			x->a[j][0] += dx;
			J->a[i][j] = (f[i](x) + minusfx->a[i][0])/dx;
			x->a[j][0] -= dx;
		}

	return J;

}



// Solves f(x) = 0 with newtons method, where f is a function
// from R^N -> R^N
// f is an array with the N coordinate functions of f,
// with x an initial guess. x should be a column vector of length N,
// and dx should be a size that's small compared to the scale of the problem.
// acc is the desired accuracy.

// Returns the number of iterations used, with the root returned in x.
// // If this count is greater than 5000, usually something bad has happened
// // and the method has failed to converge.
//

int newton(size_t n, double (*f[n])(matrix *x), matrix *x, double dx, double acc) {

	matrix *J; // Jacobian matrix	
	matrix *Dx; // Increment at each step
	matrix *z; // Test x
	matrix *minusfx; // Column vector containing -f(x)
	matrix *minusfx2; // Column vector containing -f(x) while backtracking
	matrix *R;
	double lambda;
	size_t i;
	int iterations;

	if (x->n != n || x->m != 1) {
		fprintf(stderr,"Wrong coordinate dimensions.\n");
		exit(-1);
	}

	dx = (dx <= 0) ? 1e-3 : dx;
	acc = (acc <= 0) ? 1e-6 : acc;

	J = matrix_create(n,n); // Some initializations
	R = matrix_create(n,n);
	Dx = matrix_create(n,1);
	z = matrix_create(n,1);
	minusfx = matrix_create(n,1);
	minusfx2 = matrix_create(n,1);

	for (i = 0; i < n; i++)
		minusfx->a[i][0] = -f[i](x);

	iterations = 0;
	do { // Main loop
	iterations++;
	if (unlikely(iterations > 5000)) {
		fprintf(stderr,"Warning: newton method failed to converge.\n");
		break;
	}
	// Find the next step direction Dx such that JDx = -f(x)
	// This involves first calculating the jacobian:
	J = square_jacobian(n, f, J, x, minusfx, dx);

	// Now solve J * Dx = -f(x) by QR decomposing J = QR
	// and solving R * Dx = Q^T(-f(x))
	
	QR_decomposition_prealloc(J, R); // If J is singular, we
					 // ought to do something intelligent
					 // here
	solve_linear_system(J,R,minusfx,Dx);

	// Now do backtracking line search to approximately minimize
	// ||f(x+lambda*Dx)||
	
	matrix_scale(Dx,2);
	lambda = 2.0;
	do {
		lambda /= 2.0;
		// z = x + lambda * Dx : (lambda is divided by 2 in each step)
		matrix_copy_preallocated(z, x);
		matrix_add(z,matrix_scale(Dx,0.5));
	
		for (i = 0; i < n; i++)
			minusfx2->a[i][0] = -f[i](z);
	} while (matrix_vector_norm(minusfx2) > (1-lambda/2)  *  \
		 matrix_vector_norm(minusfx) && lambda > 1.0/128.0);

	matrix_copy_preallocated(minusfx,minusfx2);
	matrix_copy_preallocated(x,z);

	} while (matrix_vector_norm(minusfx) > acc); // End main loop.

	matrix_free(J);
	matrix_free(R);
	matrix_free(Dx);
	matrix_free(z);
	matrix_free(minusfx);
	matrix_free(minusfx2);

	return iterations;

}


