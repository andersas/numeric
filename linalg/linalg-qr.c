#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "pp.h"
#include "matrix.h"
#include "linalg.h"


struct QR_worker_arg {

	size_t curr_idx;
	matrix *Q, *R;
	int command;
	size_t start, stop;
	barrier_t *b;

};

void * QR_worker(void *ptr) {

	struct QR_worker_arg * arg = (struct QR_worker_arg *) ptr;
	size_t i,j,k;
	size_t start, stop;

	while (1) {
	// Wait for commander
	barrier(arg->b);

	if (arg->command == 1) {
	start = arg->start;
	stop = arg->stop;
	i = arg->curr_idx;	
	for (j = start; j < stop; j++) {
		arg->R->a[i][j] = matrix_dot_cols(arg->Q,arg->Q,i,j);
		for (k = 0; k < arg->Q->n; k++) // Loop along rows
			arg->Q->a[k][j] -= arg->Q->a[k][i] * arg->R->a[i][j];
	}

	// Signal that we're done.
	barrier(arg->b);
	}
	
	else { // Command != 1
		pthread_exit(NULL);
	}

	}
}


// Take the QR decomposition of A in-place,
// returning Q via A. R is returned.
// A is an nxm, n>=m, matrix and R is an mxm triangular matrix.

matrix * QR_decomposition(matrix *A) {

	size_t n,m;

	n = A->n;
	m = A->m;
	if (n < m) {
		fprintf(stderr, "Cannot create orthogonal matrix from %zux%zu matrix A.\n", n, m);
		exit(-1);
	}

	return QR_decomposition_prealloc(A, matrix_create_zero(m,m));
}

// Note: lower half of R is not touched, and should be given as zero
matrix * QR_decomposition_prealloc(matrix *A, matrix *R) {

	size_t n,m,i,j,k;
	size_t n_threads = numCPUs()*2;
	pthread_t thread[n_threads];
	struct QR_worker_arg t_arg[n_threads];
	barrier_t sync_barrier;
	pthread_attr_t attr;

	n = A->n;
	m = A->m;
	if (n < m ) {
		fprintf(stderr, "Cannot create orthogonal matrix from %zux%zu matrix A.\n", n, m);
		exit(-1);
	} else if (R->n != m || R->m != m) {
		fprintf(stderr, "QR decomposition: R matrix has wrong shape\n");
		exit(-1);
	}

	// Use modified Gram-Schmidt orthogonalization

	// Do we need to parallelize?

	if (m < 4 * n_threads || numCPUs() == 1) { // no
	
	for (i = 0; i < m; i++) { // Loop over columns
		// Normalize the i'th column:
		R->a[i][i] = sqrt(matrix_dot_cols(A,A,i,i));
		if (fabs(R->a[i][i]) < 1.0e-12) {
			fprintf(stderr, "QR decomposition hit a singular matrix.\n");
			exit(-1);
		}
		for (k = 0; k < n; k++)
			A->a[k][i] = A->a[k][i]/R->a[i][i];

		// Then make all remaining columns orthogonal to column i.
		for (j = i+1; j<m;j++) {
			R->a[i][j] = matrix_dot_cols(A,A,i,j);
			for (k = 0; k < n; k++) // Loop along rows
				A->a[k][j] -= A->a[k][i] * R->a[i][j];
		}
	}
	} else { // Parallelize
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
		barrier_init(&sync_barrier, n_threads + 1);

		for (i = 0; i < n_threads; i++) {
		   t_arg[i].Q = A;
		   t_arg[i].R = R;
		   t_arg[i].command = 1;
		   t_arg[i].b = &sync_barrier;
		   pthread_create(&thread[i], &attr, QR_worker, &t_arg[i]);
		}
		pthread_attr_destroy(&attr);

		/* We now have a number of threads running, awaiting
 		 * commands. Command 0 is for thread exit,
 		 * command 1 is for running orthogonalization */

	for (i = 0; i < m; i++) { // Loop over columns
		// Normalize the i'th column:
		R->a[i][i] = sqrt(matrix_dot_cols(A,A,i,i));
		if (fabs(R->a[i][i]) < 1.0e-12) {
			fprintf(stderr, "QR decomposition hit a singular matrix.\n");
			exit(-1);
		}
		for (k = 0; k < n; k++)
			A->a[k][i] = A->a[k][i]/R->a[i][i];

		// Then make all remaining columns orthogonal to column i,
		// parallelly

		if (m - (i+1) > 4*n_threads) {	
		for (k = 0; k < n_threads; k++) {
			t_arg[k].curr_idx = i;
			t_arg[k].start = i+1 + (k*(m - (i+1)))/n_threads;
			t_arg[k].stop = i+1 + ((k+1)*(m - (i+1)))/n_threads;
		}

		// Start the orthogonalization
		barrier(&sync_barrier);

		// Then wait for it to be finished and run the
		// next column.

		barrier(&sync_barrier);

		} else { // Don't parallelize the last bit
		for (j = i+1; j<m;j++) {
			R->a[i][j] = matrix_dot_cols(A,A,i,j);
			for (k = 0; k < n; k++) // Loop along rows
				A->a[k][j] -= A->a[k][i] * R->a[i][j];
		}	
		}

		}
		for (i = 0; i < n_threads; i++) {
		   t_arg[i].command = 0;
		}
		// Ready to send command 0 to threads.
		barrier(&sync_barrier);
		barrier_destroy(&sync_barrier);

	}


	return R;
}


double QR_absdet(matrix *Q, matrix *R) {
	size_t i,m,n;

	double prod = 1.0;

	n = Q->n;
	m = Q->m;

	if (n != m) {
		fprintf(stderr, "Cannot take the determinant of a non-square matrix\n");
		exit(-1);
	} else if (R->n != m || R->m != m) {
		fprintf(stderr, "QR_absdet: R matrix has wrong dimensions).\n");
	}

	for (i = 0; i<n; i++) {
		prod *= R->a[i][i];
	}

	return prod;
}



