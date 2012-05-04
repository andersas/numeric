#!/usr/bin/env python

from math import *;

def f(x):
	if (x > 0.0):
		return 1.0;
	else:
		return 0.0;

for i in range(-10,11):
	print(i, "\t", f(1.0 * i));

