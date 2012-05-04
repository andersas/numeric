#include "data.h"
#ifndef QSPLINE
#define QSPLINE 



struct qsplines_workspace {

	struct point *data;
	double *a;
	double (*eval)(struct qsplines_workspace *state, double z);
	int count;
	int last_index;
};


struct qsplines_workspace * qspline_init(int npoints, const struct point points[]);


#endif
