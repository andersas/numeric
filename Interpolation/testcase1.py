#!/usr/bin/env python

from math import *;

def f(x):
	if (abs(x) < 0.00001):
		return 1.0;
	else:
		return sin(x)/x;

for i in range(-10,11):
	print(i, "\t", f(1.0 * i));

