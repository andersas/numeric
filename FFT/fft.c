#include <complex.h>
#include <math.h>
#include "fft.h"


// dft. Set sign to 1 (instead of -1) for inverse transform.

void dft(int N, complex *x, complex *c, int sign) {

	complex w = cexp(sign * 2 * M_PI * I / N);
	int k,n;


	for (k = 0; k < N; k++) {
		c[k] = 0;
		for (n = 0; n < N; n++) {
			c[k] += x[n]*cpow(w, n*k);
		}
		c[k] /= sqrt(N);
	}
}

void fft(int N, complex *x, complex *c, int sign) {

	int M,k,m;
	complex w;

//	if (N % 2 == 0 && N > 20) {
	if (N % 2 == 0) {
		w = cexp(sign * 2 * M_PI * I / N);

		// Split the arrays in two and perform fourier transformation
		// on the subarrays

		M = N/2;
		complex xodd[M], xeven[M], codd[M],ceven[M];

		for (m = 0; m < M; m++)
			xodd[m] = x[2*m+1],
			xeven[m] = x[2*m];
			

		fft(M, xodd, codd, sign);
		fft(M, xeven, ceven, sign);

		// Now apply twiddle factor:
		// (the sqrt 2 factor is to multiply the normalization factor
		// in the two lower transforms).
		
		for (k = 0; k < M; k++)
			c[k] = (ceven[k] + codd[k]*cpow(w,k))/M_SQRT2;
		for (k = M; k < N; k++)
			c[k] = (ceven[k-M] - codd[k-M]*cpow(w,k-M))/M_SQRT2;

	} else {
		dft(N, x, c, sign);
	}
}


