#include <stdio.h>
#include <math.h>
#include "../linalg/linalg.h"
#include "newton.h"
#include "simplex.h"

double banana(double *p) {
	double x, y;
	x = p[0];
	y = p[1];
	return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
//	return x*x+y*y;
}

double Himmelblau(double *p) {
	double x, y;
	x = p[0];
	y = p[1];
	return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}

// Mexican hat potential
double sombrero(double *p) {
	double x, y;
	x = p[0];
	y = p[1];
	
	return -10*(x*x+y*y) + (x*x+y*y)*(x*x+y*y);
}

// (1-x)=0, 10(y-x²)=0

double f0(matrix *x) {
	return 1 - x->a[0][0];
}

double f1(matrix *p) {
	double x,y;
	x = p->a[0][0];
	y = p->a[1][0];
	return 10*(y - x*x);
}

double poly(matrix *p) {
	return p->a[0][0]*p->a[0][0] - 101;
}

double circ1(matrix *p) {
	double x,y,z;
	x = p->a[0][0];
	y = p->a[1][0];
	z = p->a[2][0];
	return (x-5)*(x-5) + y*y + z*z - 25;
}
double circ2(matrix *p) {
	double x,y,z;
	x = p->a[0][0];
	y = p->a[1][0];
	z = p->a[2][0];
	return x*x + y*y + z*z - 9;
}
double circ3(matrix *p) {
	double x,y,z;
	x = p->a[0][0];
	y = p->a[1][0];
	z = p->a[2][0];
	return (x-1)*(x-1) + (y+2)*(y+2) + (z+4)*(z+4) - 25;
}
int main() {
	double (*first[2])(matrix *x) = {&f0, &f1};
	double (*sqrtsolve[1])(matrix *x) = {&poly};
	double (*circsyst[3])(matrix *x) = {&circ1, &circ2, &circ3};
	matrix *init = matrix_create(2,1);
	double f;
	int c;
	init->a[0][0] = 2;
	init->a[1][0] = 0;

	printf("- - - - - - - - - - - - - -\n");
	printf("           Start           \n");
	printf("- - - - - - - - - - - - - -\n");
	printf("Finding minimum of banana function...\n");
	f = downhill_simplex_method(2,&banana, init, 1, 1e-6);
	printf("banana(%f, %f) = %f.\n", init->a[0][0], init->a[1][0], f);

	init->a[0][0] = -50;
	init->a[1][0] = -50;
	printf("Finding minimum of Himmelblau function with starting a simplex that's a square centered on 0...\n");
	f = downhill_simplex_method(2,&Himmelblau, init, 100, 1e-6);
	printf("Himmelblau(%f, %f) = %f.\n", init->a[0][0], init->a[1][0], f);
	printf("(there are several, but the function rises extremely fast away from them)\n");

	printf("Seeing what happens with the mexican hat potential function...\n");
	init->a[0][0] = -5;
	init->a[1][0] = -5;
	f = downhill_simplex_method(2,&sombrero, init, 10, 1e-6);
	printf("sombrero(%f, %f) = %f.\n", init->a[0][0], init->a[1][0], f);
	printf("x^2+y^2 = %f, should be 5.\n",init->a[0][0]*init->a[0][0]+init->a[1][0]*init->a[1][0]);
	printf("Behold! Broken symmetry!\n");

	printf("\n- - - - - - - - - - - - - -\n");
	printf("Switching to Newtons method..\n");
	printf("- - - - - - - - - - - - - -\n\n");

	printf("Solving (1-x)=0, 10(y-x²)=0...\n");
	
	init->a[0][0] = 5;
	init->a[1][0] = 5;
	c = newton(2, first, init,1e-3, 1e-6);
	printf("x = %f, y = %f solves the system (after %i iterations).\n",init->a[0][0],init->a[1][0], c);

	printf("Calculating the square root of 101 by solving\n");
	printf("x^2 - 101 = 0\n");

	matrix_free(init);
	init = matrix_create(1,1);
	init->a[0][0] = 100;
	c = newton(1,sqrtsolve, init, 1e-3, 1e-6);
	printf("sqrt(101) = %f (in fact %f) after %i iterations.\n",init->a[0][0],sqrt(101), c);

	printf("\nI wonder how GPS works...\n");
	printf("Consider 3 spheres given by centrum and radius:\n");
	printf("(5,0,0), radius 5,\n");
	printf("(0,0,0), radius 3,\n");
	printf("(1,-2,-4), radius 5.\n");
	printf("Finding intersection by solving the system\n");
	printf("(x-5)^2 + y^2 + z^2 - 5^2 = 0,\n");
	printf("x^2 + y^2 + z^2 - 3^2 = 0,\n");
	printf("(x-1)^2 + (y+2)^2 + (z+4)^2 - 5^2 = 0...\n");

	matrix_free(init);
	init = matrix_create(3,1);
	init->a[0][0] = 0;
	init->a[1][0] = 0;
	init->a[2][0] = 0;
	c = newton(3,circsyst, init, 1e-3, 1e-6);
	printf("Intersection at (%f, %f, %f) found after %i iterations.\n", init->a[0][0], init->a[1][0], init->a[2][0], c);

return 0;
}
