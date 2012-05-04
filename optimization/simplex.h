#ifndef _SIMPLEX
#define _SIMPLEX
#include <stdio.h>
#include "../linalg/linalg.h"

typedef struct simplex {
	double (*f)(double *p); // Pointer to function we creep on
	double *fp; // array of function values in each point
	matrix *centroid; // center of gravity for all points but Pmax_idx
	matrix *points; // Matrix with each point as a row vector.
	size_t dim; // Dimension of the simplex

} simplex;

// The downhill simplex method function.
// dim is the dimension of the function, f is the function
// and init is a starting point for the simplex (column or row vector).
// The simplex is a dim dimensional cube with
// one corner being init and side lengths range.
// tol is the tolerance. After the simplex has shrunk below
// this size, the lowest point is returned via init.

// Returns the minimum value found and updates init to be the
// point where this value is found.

double downhill_simplex_method(size_t dim, double (*f)(double *p), matrix *init, double range, double tol);

// Create a simplex centered on init with other points at most range from init.
simplex *simplex_create(size_t dim, double (*f)(double *p), matrix *init, double range);
void simplex_free(simplex *s);

void simplex_sort_points(simplex *s);
void simplex_update_highest(simplex *s);
void simplex_update_centroid(simplex *s);


// Operations on simplexes: reflection, expansion, contraction and reduction.
// The first three are conditional operations that must be accepted
// before being applied to the actual simplex.

// Pre is a 1xdim matrix that will contain the reflection.
// f(Pre) is returned.
double simplex_get_reflection(simplex *s, matrix *Pre);
double simplex_get_expansion(simplex *s, matrix *Pex);
int simplex_get_contraction(simplex *s, matrix *Pco);

// Update the simplex
void simplex_accept(simplex *s, matrix *Phi_new, double fPhi_new);
void simplex_reduce(simplex *s);
double simplex_size(simplex *s);

#endif
