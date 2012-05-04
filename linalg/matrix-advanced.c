#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "pp.h"

#define min(a,b) ((a < b) ? a : b)
#define max(a,b) ((a > b) ? a : b)

matrix *matrix_transpose(matrix *mat) {

	size_t i,j,n,m;
	matrix *buf;

	n = mat->n;	
	m = mat->m;	

	// This is a temporary stub that just works by
	// doing an out-of-place transposition
	// in the naive way. For better effeciency, implement
	// one of the in-place algorithms suggested on wikipedia.
	// For now this will have to suffice.


	buf = matrix_create(m, n);

	for (i = 0; i < n; i++)
		for(j = 0; j < m; j++)
			buf->a[j][i] = mat->a[i][j];

	free(mat->root);
	free(mat->a);

	mat->a = buf->a;
	mat->root = buf->root;

	mat->n = buf->n;
	mat->m = buf->m;

	free(buf);
	return mat;
}

struct matrix_multiply_worker_arg {

	matrix *C, *A,*B;

	int work_by_rows;
	size_t start, stop;	

};

void * matrix_multiply_worker(void * arg) {

	struct matrix_multiply_worker_arg *p = (struct matrix_multiply_worker_arg *) arg;
	matrix *C, *A, *B;
	size_t i,j;

	C = p->C;
	A = p->A;
	B = p->B;

	if (p->work_by_rows) {
		for (i = p->start; i < p->stop; i++)
		for (j = 0; j < C->m; j++)
			C->a[i][j] = matrix_dot_product(A,i,B,j);
	} else {
		for (i = 0; i < C->n; i++)
		for (j = p->start; j < p->stop; j++)
			C->a[i][j] = matrix_dot_product(A,i,B,j);
	}

	return NULL;
}



matrix *matrix_multiply(matrix *A, matrix *B) {

	if (A->m != B->n) {
		fprintf(stderr, "Attempting to multiply a %zux%zu matrix by a %zux%zu matrix. Bailing!\n", A->n, A->m, B->n, B->m);
		exit(-1);
	}

	
	return matrix_multiply_prealloc(matrix_create(A->n, B->m),A,B);

}

matrix *matrix_multiply_prealloc(matrix *C, matrix *A, matrix *B) {

	// Implement the naive matrix multiplication algorithm
	// in parallel. If the time allows it,
	// In future try to implement the Strassen Matrix Multiplication algorithm
	// that runs O(n^2.8) insted. Strassen is also parallelizable.
	
	size_t n_threads = numCPUs()*2;
	struct matrix_multiply_worker_arg thread_arg[n_threads];
	size_t i,j;
	pthread_t thread[n_threads];
	pthread_attr_t attr;

	if (A->m != B->n) {
		fprintf(stderr, "Attempting to multiply a %zux%zu matrix by a %zux%zu matrix. Bailing!\n", A->n, A->m, B->n, B->m);
		exit(-1);
	} else if (C->n != A->n || C->m != B->m) {
		fprintf(stderr, "Wrong dimensions on preallocated result matrix.\n");
		exit(-1);
	}

	if (A->n * A->m * B->m < 100000 || max(C->n, C->m) < numCPUs() || numCPUs() == 1) {
		// Do less than 10^5 double multiplications.
		// It's probably faster to just do it than
		// spawn new threads.

		for (i = 0; i < C->n; i++)
		for (j = 0; j < C->m; j++)
			C->a[i][j] = matrix_dot_product(A,i,B,j);

	} else {
		// Let's spread the load over available cores.
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		for (i = 0; i < n_threads; i++) {

		thread_arg[i].C = C;
		thread_arg[i].A = A;
		thread_arg[i].B = B;

		if (C->n >= C->m) {
			// Loop over rows
			thread_arg[i].work_by_rows = 1;
			thread_arg[i].start = i*(C->n)/(n_threads);
			thread_arg[i].stop = min((i+1)*C->n/(n_threads),C->n);
		} else {
			// Loop over columns
			thread_arg[i].work_by_rows = 0;
			thread_arg[i].start = i*(C->m)/(n_threads);
			thread_arg[i].stop = min((i+1)*C->m/(n_threads),C->m);
		}

		// Dispatch a worker.
		pthread_create(&thread[i], &attr, &matrix_multiply_worker, (void *) &thread_arg[i]);

		}
	pthread_attr_destroy(&attr);

		for (i = 0; i < n_threads; i++) {
			pthread_join(thread[i], NULL);
		}

	}
	
	return C;
}

