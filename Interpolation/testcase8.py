#!/usr/bin/env python

from math import *;
from random import *;

def f(x):
	if (x > 0.0):
		return 1.0 + gauss(0,0.02);
	else:
		return 0.0 + gauss(0,0.02);

for i in range(-10,11):
	print(i, "\t", f(1.0 * i));

