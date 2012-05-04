#include "step.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Add a constant times a vector to a vector
inline void addscaled(int len, double dest[len], double v, const double src1[len], double w, const double src2[len]) {
	int i;
	for (i = 0; i < len; i++)
		dest[i] = v*src1[i] + w*src2[i];
}

// returns error estimate, values in y0
double rkstep1(void (*f)(double x, const double y[], double result[]), int dim, double x0, double y0[dim], double h) {

	// Use Eulers method to calculate f and use two half steps
	// to estimate the error. 

	double k0[dim];  // Deriviate at x0
	double k12[dim]; // Deriviate at half step
	double y1[dim]; // Increment in F at full step
	//double y0[dim];  // Increment in F at two half steps
	double error;
	int i;

	f(x0, y0, k0);
	addscaled(dim, y1, 1.0, y0, h, k0);
	f(x0 + h/2.0, y1, k12);
	addscaled(dim, k0, h/2.0, k0, h/2.0, k12);
	addscaled(dim, y0, 1.0, y0, 1.0, k0);

	error = 0.0;
	for (i = 0; i < dim; i++) {
		error += (y1[i]-y0[i])*(y1[i]-y0[i]);
	}
	error = sqrt(error);

	return error;

}



// Use a fourth order method. 
void rksubstep4(void (*f)(double x, const double y[], double result[]), int dim, double x0, double y0[dim], double h, const double initial[dim]) {


	double Kimm[dim]; // K intermediate result.
	double yimm[dim];
	double k[dim]; // finished K
	int i;

		// k = Kimm = initial
		addscaled(dim, k, 1.0, initial, 0.0, initial);
		addscaled(dim, Kimm, 1.0, initial, 0.0, initial);
		// yimm = y0 + h/2 * Kimm
		addscaled(dim, yimm, 1.0, y0, h/2.0, Kimm);
		// k += 2 * (Kimm = f(x0 + h/2.0, yimm));
		f(x0 + h/2.0, yimm, Kimm);
		addscaled(dim, k, 1.0, k, 2.0, Kimm);
		// k += 2 * (Kimm = f(x0 + h/2.0, yimm));
		f(x0 + h/2.0, yimm, Kimm);
		addscaled(dim, k, 1.0, k, 2.0, Kimm);
		// k[i] += (Kimm[i] = f(x0 + h, yimm));
		f(x0 + h, yimm, Kimm);
		addscaled(dim, k, 1.0, k, 1.0, Kimm); 
		//k /= 6.0;
		for (i = 0; i < dim; k[i++] /= 6.0);
		//y0 += h*k;
		addscaled(dim,y0,1.0, y0, h, k);

}

// This is an intermediate driver that takes two half steps and
// one full step and compare to get the error.
double rkstep4(void (*f)(double x, const double y[], double result[]), int dim, double x0, double y0[dim], double h) {

	double ytmp[dim];
	double initial[dim];
	double error;
	int i;

	memcpy(ytmp, y0, dim*sizeof(double));

	f(x0,y0,initial);

	rksubstep4(f,dim,x0,ytmp,h,initial); // full step (update ytmp)
	// Then two half steps:
	rksubstep4(f,dim,x0,y0,h/2.0,initial); // Substep changes y0
	f(x0+h/2.0, y0, initial);
	rksubstep4(f,dim,x0,y0,h/2.0,initial); // (and so we also change y0)

	error = 0.0;
	for (i = 0; i < dim; i++) {
		error += (ytmp[i]-y0[i])*(ytmp[i]-y0[i]);
	}

	return sqrt(error);
}


