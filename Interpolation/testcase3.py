#!/usr/bin/env python

from math import *;

def f(x):
	return exp(-1.0/20*x**2);

for i in range(-10,11):
	print(i, "\t", f(1.0 * i));

