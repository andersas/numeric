#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import *;
from random import *;

print("#", 10, ",\t", 3);
for t in range(0,10):
	print(15.2*(cos(pi/5.0*t) + gauss(0,0.06))+2.4, "\t", 15.2*(sin(pi/5.0*t) + gauss(0,0.06))-3.6, "\t", sqrt(2*0.06**2), "\n");


