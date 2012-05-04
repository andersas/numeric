#!/usr/bin/env python

from math import *;
from random import *;

def f(x):
	if (abs(x) < 0.00001):
		ret = 1.0;
	else:
		ret = sin(x)/x;
	return ret + gauss(0,0.02);

for i in range(-10,11):
	print(i, "\t", f(1.0 * i));

