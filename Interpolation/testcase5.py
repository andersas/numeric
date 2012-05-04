#!/usr/bin/env python

from math import *;

def f(x):
	return 1/abs(x-2.5);

for i in range(-10,11):
	print(i, "\t", f(1.0 * i));

