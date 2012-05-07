#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "matrix.h"

double matrix_dot_product_part(double **matrix, size_t col, double *row, size_t a, size_t b) {
	double sum = 0.0;
	size_t i;

	// Pairwise summation algorithm for dot product
	// between column (non-contingous) and a row (contingous) posize_ted to by 
	// row.

	if ((b - a) < 300) {
		for (i = a; i<b; i++) {
			sum += row[i] * matrix[i][col];
		}
	} else {
		sum += matrix_dot_product_part(matrix,col,row,a,(b+a)/2);
		sum += matrix_dot_product_part(matrix,col,row,(b+a)/2,b);
	}
	return sum;
}

// Wrapper for the above function:
double matrix_dot_product(matrix *A, size_t row, matrix *B, size_t col) {
	if (A->m != B->n) {
		fprintf(stderr, "Invalid dot product!\n");
		exit(-1);
	}
	return matrix_dot_product_part(B->a, col, A->a[row], 0,A->m);
}

// Take the dot product of two vectors using pairwise summation

double matrix_single_dot_part(double *v, double *w, size_t a, size_t b) {
	double sum = 0.0;
	size_t i;

	if ((b - a) < 300) {
		for (i = a; i<b; i++) {
			sum += v[i]*w[i];
		}
	} else {
		sum += matrix_single_dot_part(v,w,a,(b+a)/2);
		sum += matrix_single_dot_part(v,w,(b+a)/2,b);
	}
	return sum;
}

// Wrapper
double matrix_single_dot_product(size_t len, double *v, double *w) {
	return matrix_single_dot_part(v, w, 0, len);
}


double matrix_dot_cols_part(double **matrix1, double **matrix2, size_t col1, size_t col2, size_t a, size_t b) {
	double sum = 0.0;
	size_t i;

	// Pairwise summation algorithm for dot product
	// between two columns in matrices (non-contingous)

	if ((b - a) < 300) {
		for (i = a; i<b; i++) {
			sum += matrix1[i][col1] * matrix2[i][col2];
		}
	} else {
	   sum += matrix_dot_cols_part(matrix1,matrix2,col1,col2,a,(b+a)/2);
	   sum += matrix_dot_cols_part(matrix1,matrix2,col1,col2,(b+a)/2,b);
	}
	return sum;
}

// Wrapper
double matrix_dot_cols(matrix *A, matrix *B, size_t colA, size_t colB) {
	if (A->n != B->n) {
		fprintf(stderr, "Cowardly refusing to take the dot product of vectors of different length.\n");
		exit(-1);
	}

	return matrix_dot_cols_part(A->a, B->a, colA, colB, 0, A->n);
}

double matrix_vector_norm(matrix *v) {

	if (v->n == 1) {
		// row vector
		return sqrt(matrix_single_dot_product(v->m,v->root,v->root));
	} else if (v->m == 1) {
		// column vector
		return sqrt(matrix_single_dot_product(v->n,v->root,v->root));
	} else {
		// Uh-oh
		if (v->n == v->m && v->n == 1) // both
			return fabs(v->a[0][0]);

		fprintf(stderr,"Refusing to take norm of non-vector.\n");
		exit(-1);
	}
}

// Allocate a new matrix
matrix *matrix_create(size_t n, size_t m) {

	size_t i;
	matrix *mat = (matrix *) malloc (sizeof(matrix));
	double **q = (double **) malloc(n*sizeof(double));

	mat->n = n;
	mat->m = m;
	mat->a = q;
	mat->root = (double *) malloc(n*m*sizeof(double));

	for (i = 0; i < n; i++) {
		q[i] = &(mat->root[m*i]);
	}

	return mat;
}

// Allocate a new matrix filled with zeros (using calloc instead of malloc)
matrix *matrix_create_zero(size_t n, size_t m) {

	size_t i;
	matrix *mat = (matrix *) malloc (sizeof(matrix));
	double **q = (double **) malloc(n*sizeof(double));

	mat->n = n;
	mat->m = m;
	mat->a = q;
	// This assumes IEEE 754 double precision numbers
	mat->root = (double *) calloc(n*m, sizeof(double));

	for (i = 0; i < n; i++) {
		q[i] = &(mat->root[m*i]);
	}

	return mat;
}

// Allocate a matrix and fill it with random numbers
// from GSL mersenne twister
matrix *matrix_create_random(size_t n, size_t m) {
	size_t i,j;
	unsigned long int seed;
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
	matrix *mat = matrix_create(n,m);

	FILE *F;
	if ((F = fopen("/dev/urandom", "r")) == NULL) {
		perror("fopen");
		exit(-1);
	}
	if (!fread(&seed, sizeof(unsigned long int), 1, F)) {
		perror("fread");
		exit(-1);
	}
	fclose(F);
	gsl_rng_set(r,seed);

	for (i=0; i<n; i++)
		for (j=0; j<m; j++)
			mat->a[i][j] = gsl_rng_uniform(r);

	gsl_rng_free(r);

	return mat;
	}

// Add two matrices (store result in A)
matrix *matrix_add(matrix *A, matrix *B) { // A = A + B

	size_t i,j;
	if (A->n != B->n || A->m != B->m) {
		fprintf(stderr, "Trying to add a matrix of dimension %zux%zu with a matrix of dimension %zux%zu\nCowardly aborting!\n", A->n,A->m,B->n,B->m);
		exit(-1);
	}

	for (i = 0; i < A->n; i++)
		for (j = 0; j < A->m; j++)
			A->a[i][j] += B->a[i][j];

	return A;
}

// Subtract two matrices (store result in A)
matrix *matrix_sub(matrix *A, matrix *B) { // A = A - B

	size_t i,j;
	if (A->n != B->n || A->m != B->m) {
		fprintf(stderr, "Trying to subtract a matrix of dimension %zux%zu with a matrix of dimension %zux%zu\nCowardly aborting!\n", A->n,A->m,B->n,B->m);
		exit(-1);
	}

// Rows may be switched, so we don't use root
	for (i = 0; i < A->n; i++)
		for (j = 0; j < A->m; j++)
			A->a[i][j] -= B->a[i][j];

	return A;
}


// Multiply a matrix by a scalar (store result in A)
matrix *matrix_scale(matrix *A, double b) {
	size_t i;

	size_t stop = A->n * A->m;

	for (i = 0; i < stop; i++)
		A->root[i] *= b;

/*	for (i = 0; i < A->n; i++)
		for (j = 0; j < A->m; j++)
			A->a[i][j] *= b;
*/
	return A;
}

// Scale single column
matrix *matrix_scale_col(matrix *mat, size_t col, double b) {
	size_t i;

	if (mat->m < col) {
		fprintf(stderr, "Refusing to scale column %zu of and %zux%zu matrix!\n", col, mat->n, mat->m);
		exit(-1);
	}

	for (i = 0; i < mat->n; i++)	
		mat->a[i][col] *= b;

	return mat;
}

// same for row
matrix *matrix_scale_row(matrix *mat, size_t row, double b) {
	size_t i;

	if (mat->n < row) {
		fprintf(stderr, "Refusing to scale row %zu of and %zux%zu matrix!\n", row, mat->n, mat->m);
		exit(-1);
	}

	for (i = 0; i < mat->m; i++)	
		mat->a[row][i] *= b;

	return mat;
}


// Make a diagonal matrix with the diagonal elements in a vector
matrix *matrix_create_diagonal(size_t n, matrix *diagonal) {
	matrix *mat = matrix_create_zero(n,n);
	return matrix_set_diagonal(mat, diagonal);
}

matrix *matrix_create_diagonal_scalar(size_t n, double diagonal) {
	return matrix_set_diagonal_scalar(matrix_create_zero(n,n),diagonal);
}

matrix *matrix_identity(size_t n) {
	return matrix_create_diagonal_scalar(n, 1.0);
}

// Copy diagonal elements of square matrix to storage area posize_ted to by storage
matrix *matrix_get_diagonal(matrix *mat, matrix *storage) {
	size_t i;
	if (mat->n != mat->m) {
		fprintf(stderr,"Trying to get diagonal of non-square matrix -- Aborting!\n");
		exit(-1);
	} if (storage->n != mat->n || storage->m != 1) {
		fprintf(stderr,"Invalid size for storage matrix.\n");
		exit(-1);
	}
	for (i = 0; i < mat->n; i++) {
		storage->a[i][0] = mat->a[i][i];
	}

	return storage;	
}

// The opposite of above
matrix *matrix_set_diagonal(matrix *mat, const matrix *storage) {
	size_t i;
	if (mat->n != mat->m) {
		fprintf(stderr,"Trying to set diagonal of non-square matrix -- Aborting!\n");
		exit(-1);
	} if (storage->n != mat->n || storage->m != 1) {
		fprintf(stderr,"Invalid size for storage matrix.\n");
		exit(-1);
	}

	for (i = 0; i < mat->n; i++) {
		mat->a[i][i] = storage->a[i][0];
	}
	return mat;
}

matrix *matrix_set_diagonal_scalar(matrix *mat, double diagonal) {
	size_t i;
	if (mat->n != mat->m) {
		fprintf(stderr,"Trying to set diagonal of non-square matrix -- Aborting!\n");
		exit(-1);
	}

	for (i=0; i<mat->n; i++)
		mat->a[i][i] = diagonal;

	return mat;
}

matrix *matrix_zero(matrix *mat) {

	// This assumes IEEE 754 double precision numbers
	memset(mat->root, 0, (mat->n) * (mat->m) * sizeof(double));
	return mat;
}

// Notice how this function breaks assumptions on how
// A->root works ! ! !
matrix *matrix_swap_rows(matrix *A, size_t i, size_t j) {
	double *tmp;

	tmp = A->a[i];
	A->a[i] = A->a[j];
	A->a[j] = tmp;
	return A;
}


matrix *matrix_swap_cols(matrix *A, size_t i, size_t j) {
	double tmp;

	size_t k;

	for (k = 0; k < A->n; k++) {
		tmp = A->a[k][i];
		A->a[k][i] = A->a[k][j];
		A->a[k][j] = tmp;
	}

	return A;
}


// Call this to dispose of used matrices
void matrix_free(matrix * mat) {
	free(mat->root);
	free(mat->a);
	free(mat);
}

matrix *matrix_column_from_array(size_t len, const double v[]) {
	size_t i;
	matrix *mat = matrix_create(len,1);

	for (i=0; i < len ; i++)
		mat->a[i][0] = v[i];

	return mat;
}

matrix *matrix_row_from_array(size_t len, const double v[]) {
	size_t i;
	matrix *mat = matrix_create(1,len);

	for (i=0; i < len ; i++)
		mat->a[0][i] = v[i];

	return mat;
}


matrix *matrix_copy(matrix * mat) {
	return matrix_copy_preallocated(matrix_create(mat->n, mat->m), mat);
}
matrix *matrix_clone(matrix * mat) {
	return matrix_copy(mat);
}

matrix *matrix_copy_preallocated(matrix *dest, matrix *src) {
	size_t i;
	if (dest->n != src->n || dest->m != src->m) {
		fprintf(stderr, "Copying matrix to destination of incompatible dimensions.\n");
		exit(-1);
	}

	for (i = 0; i < src->n; i++)
		memcpy(dest->a[i], src->a[i], src->m*sizeof(double));

	return dest;
}
matrix *matrix_clone_preallocated(matrix *dest, matrix *src) {
	return matrix_copy_preallocated(dest, src);
}

matrix *matrix_store(matrix *mat, char *filename) {

	size_t i,j;
	FILE *f;
	if ((f = fopen(filename, "w")) != NULL) {

		fprintf(f, "# %zu, %zu\n", mat->n, mat->m);

	for (i = 0; i < mat->n; i++) {
		
		fprintf(f, "%f", mat->a[i][0]);
		for (j = 1; j < mat->m; j++) {
			fprintf(f, "\t%f", mat->a[i][j]);
		}
		fprintf(f, "\n");
	}

	fclose(f);
	return mat;

	} else {
		perror("fopen");
		exit(-1);
	}
}

matrix * matrix_load(char *filename) {

	size_t i,j;
	matrix *mat;
	FILE *f;

	if ((f = fopen(filename, "r")) != NULL) {

		if (fscanf(f, "# %zu , %zu\n", &i, &j) != 2) {
			fprintf(stderr, "Could not parse file %s.\n", filename);
			exit(-1);
		}
		mat = matrix_create(i,j);
	
		for (i = 0; i < mat->n; i++) {
			for (j = 0; j < mat->m; j++) {
				if (fscanf(f, "%lf", &(mat->a[i][j])) != 1) {
				fprintf(stderr,"Could not parse file %s.\n", filename);
				exit(-1);
				}
			}
		}
	
	fclose(f);

	return mat;
	} else {
		perror("fopen");
		exit(-1);
	}
} 

