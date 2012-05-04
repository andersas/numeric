#include <stdio.h>
#include <stdlib.h>
#include "../linalg/linalg.h"
#include "simplex.h"



double downhill_simplex_method(size_t dim, double (*f)(double *p), matrix *init, double range, double tol) {

	simplex *s = simplex_create(dim, f, init, range);
	matrix *buf = matrix_create(1,dim);
	double reflected, highest, lowest, expanded, contracted;
	size_t i, iteration_count;

	// This is the exact algorithm presented in the lecture notes
	iteration_count = 0;
	do {
	iteration_count++;
		if (iteration_count >= 5000) {
			fprintf(stderr,"Warning: downhill_simplex_method gave up after %zu iterations.\n", iteration_count);
			break;
		}
		highest = s->fp[0];
		lowest = s->fp[dim]; // there are dim + 1 entries
		reflected = simplex_get_reflection(s,buf);
		if (reflected < highest) {
			simplex_accept(s, buf, reflected);
//			simplex_update_centroid(s);
			if (reflected < lowest) {
				expanded = simplex_get_expansion(s,buf);
				if (expanded < reflected) {
					simplex_accept(s,buf,expanded);
				}
			}
		
		} else { // Reflection didn't work, try contraction instead

			contracted = simplex_get_contraction(s,buf);
			if (contracted < highest) {
				simplex_accept(s,buf,contracted);
			} else {
				simplex_reduce(s);
			}
		}
		simplex_update_centroid(s);
	} while (simplex_size(s) > tol);

	// Settled with an OK simplex. Now copy the lowest point of the simplex 
	// back into init.
	
	for (i = 0; i < dim; i++) {
		init->root[i] = s->points->a[dim][i];
	}

	matrix_free(buf);
	simplex_free(s);

	printf("Downhill simplex finished after %zu iterations.\n", iteration_count);

	return f(init->root);
} 



inline void swapd(double *a, double *b) {
	double tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

// Create a simplex centered on init with other points at most range from init.

simplex *simplex_create(size_t dim, double (*f)(double *p), matrix *init, double range) {
	size_t n,m;

	if ((init->n != 1 || init->m != dim) && (init->m != 1 || init->n != dim) ) {
		fprintf(stderr,"Wrong initial vector given to simplex_create.\n");
		exit(-1);
	}

	simplex *s = (simplex *) malloc(sizeof(simplex));
	s->fp  = (double *) malloc((dim+1)*sizeof(double));
	s->centroid = matrix_create(1,dim);
	s->points = matrix_create(dim+1,dim); // dim + 1 points, in dim dimensional space
	s->f = f;
	s->dim = dim;

	// Fill out the points
	for (n = 0; n < dim; n++)
		s->points->a[0][n] = init->root[n];

	for (n = 1; n < dim + 1; n++)
	  for (m = 0; m < dim; m++)
	    // Could be random instead of range out in one dimension
	     s->points->a[n][m] = s->points->a[0][m] + ((n-1 == m) ? range : 0.0);


	// Evaluate f in each point and sort indexes by
	// largest value (use selection sort, as we're probably
	// dealing with relatively small lists)
	
	for (n = 0; n < dim + 1; n++) {
		s->fp[n] = f(s->points->a[n]);
	}

	simplex_sort_points(s);

	// Then calculate centroid:
	
	simplex_update_centroid(s);

	return s;
}

void simplex_sort_points(simplex *s) {
	size_t n,m,dim;
	dim = s->dim;

	// selection sort is used since we're probably dealing with
	// relatively few dimensions
	for (n = 0; n < dim; n++)
		for (m = n+1; m < dim + 1; m++)
			if (s->fp[n] < s->fp[m]) {
				swapd(&(s->fp[n]), &(s->fp[m]));
				matrix_swap_rows(s->points, n,m); 
			}
}

void simplex_update_highest(simplex *s) {
	size_t i, dim;
	dim = s->dim;

	// Find the new highest point (all points are sorted except the first)
	
	for (i = 1; i < dim + 1; i++)
		if (s->fp[i] > s->fp[i-1]) {
			swapd(&(s->fp[i]), &(s->fp[i-1]));
			matrix_swap_rows(s->points, i,i-1);
		} else
			break;
}


void simplex_update_centroid(simplex *s) {
size_t n,m,dim;
dim = s->dim;
	for (n = 0; n < dim; n++) {  // Loop over coordinates
		s->centroid->a[0][n] = 0.0;
		for (m = 1; m < dim + 1; m++) { // Loop over points,
						// excluding top point
			s->centroid->a[0][n] += s->points->a[m][n];
		}
		s->centroid->a[0][n] /= dim;
	}
}

void simplex_free(simplex *s) {
	matrix_free(s->points);
	matrix_free(s->centroid);
	free(s->fp);
	free(s);
}


// Operations on simplexes: reflection, expansion, contraction and reduction.
// The first three are conditional operations that must be accepted
// before being applied to the actual simplex.


// Pre is a 1xdim matrix that will contain the reflection.
// f(Pre) is returned.
double simplex_get_reflection(simplex *s, matrix *Pre) {

	size_t i;

	// Phi -> Pre = Pce + (Pce - Phi)
	// = 2Pce - Phi
	for (i = 0; i < s->dim; i++) {
		Pre->a[0][i] = 2*s->centroid->a[0][i] - s->points->a[0][i];
	}

	return s->f(Pre->a[0]);
}
double simplex_get_expansion(simplex *s, matrix *Pex) {
	
	size_t i;

	// Phi -> Pex = Pce + 2(Pce - Phi)
	// = 3Pce - 2Phi
	for (i = 0; i < s->dim; i++)
	  Pex->a[0][i] = 3.0*s->centroid->a[0][i] - 2 * s->points->a[0][i];

        return s->f(Pex->a[0]);
}

int simplex_get_contraction(simplex *s, matrix *Pco) {

	size_t i;
	// Phi -> Pco = Pce + 1/2(Phi - Pce)
	// = (Phi + Pce)/2
	
	for (i = 0; i < s->dim; i++)
	  Pco->a[0][i] = (s->points->a[0][i] + s->centroid->a[0][i])/2.0;
        
	return s->f(Pco->a[0]);
}

// Update the simplex
void simplex_accept(simplex *s, matrix *Phi_new, double fPhi_new) {
	size_t i;
	s->fp[0] = fPhi_new;
	for (i = 0; i < s->dim; i++) // Loop over coordinates
		s->points->a[0][i] = Phi_new->a[0][i];
	simplex_update_highest(s);
//	simplex_update_centroid(s);
}

void simplex_reduce(simplex *s) {

	size_t i,j;
	size_t dim = s->dim;
	// Pk, k!=lo -> 1/2(Pk + Plo), Plo is in row dim in s->points

	for (i = 0; i < dim; i++) // Loop over points
	for (j = 0; j < dim; j++) // Loop over coordinates
	{
	 s->points->a[i][j] = (s->points->a[i][j] + s->points->a[dim][j]) / 2.0;
	}
	for (i = 0; i < dim; i++)
	 s->fp[i] = s->f(s->points->a[i]); // Reevaluate f in all points but lowest

	simplex_sort_points(s);
//	simplex_update_centroid(s);
}


double simplex_size(simplex *s) {

	size_t i,j;
	double norm;

	matrix *outer = matrix_create(1,s->dim);
	matrix *inner = matrix_create(1,s->dim);
	
	for (i = 0; i < s->dim; i++) {
		for (j = 0; j < s->dim; j++)
			inner->a[0][j] = (s->points->a[i][j] - s->points->a[s->dim][j]);
		outer->a[0][i] = matrix_vector_norm(inner);
	}

	norm = matrix_vector_norm(outer);

	matrix_free(inner);
	matrix_free(outer);

	return norm;
}


