#!/usr/bin/env python

# Make the window function
for i in range(4096):
	if i in range(1024,3072):
		print 1.0
	else:
		print 0.0

