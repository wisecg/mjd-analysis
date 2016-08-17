#!/usr/bin/env python  
# sampling.py
# Clint Wiseman, USC

import random, scipy
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.stats import gamma

# distribution to be sampled from
def f(x, p):
    mu = p[0]
    sig = p[1]
    lmb = p[2]

    # gaussian
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# integrate the function (returns a numpy list)  
def F(x,p,lo,hi):
    res = np.zeros_like(x)
    norm = 0
    for i,val in enumerate(x):
        y,err = integrate.quad(f,lo,val,args=(p))
        res[i]=y
        norm=res[i]
    return res/norm

# lookup the x value[i,0] of 'table' given the y value[i,1] (find inverse)
def lookup(y,table,pts):

	rows,cols = table.shape
	step = abs(table[1,0]-table[0,0])
	foundVal = np.array([0,0])

	for i,val in enumerate(table):
		# print "Finding y=%.2f - [%i,1](y): %.3f  [%i,0](x): %.3f  y-table[%i,1]=%.3f" % (y,i,table[i,1],i,table[i,0],i,y-table[i,1])

		if abs(table[i,1] - y) < step:
			# print "\nLooking for %.3f, found %.3f with x-val %.3f." % (y,table[i,1],table[i,0])
			# print "y: %.3f  i: %0.1f  table[i,0](x): %0.4f  abs(table[i,1]-y)= %.4f" % (y,i,table[i,0],abs(table[i,1]-y))
			foundVal = val

		if (table[i,1] - y) > step:
			# print "Final value: y = %.3f  x = %.3f" % (i,foundVal[1],foundVal[0])
			return foundVal[0]
	
	return -999 # fail to match

def main():
	
	pts = 5000
	lo = -3
	hi = 3
	p = [0,0.3,5] # mu,sig 

	print "Generating lists ..."

	linX = np.linspace(lo,hi,pts)
	rndX = np.random.rand(pts)

	func = f(linX,p)
	intg = F(linX,p,lo,hi)
	
	plt.plot(linX,func,"r-")	# function to be sampled from 
	plt.plot(linX,intg,"g-")	# integral of the function

	# combine two lists into one "pts*2"-sized matrix
	table = np.column_stack((linX,intg))

	# lookup and histogram result
	print "Looking up x-values ..."
	resultX = np.array([])
	for i,val in enumerate(rndX):
		res = lookup(val,table,pts)
		if res != -999:
			if (i%500 == 0): print "%.0f%% done." % (100*(i/float(pts)))
			resultX = np.append(resultX,[res])
	print resultX.shape
	plt.hist(resultX,50,normed=1)	# histogram of the sampling of the integral of the function.
	
	plt.show()

if __name__ == '__main__':
    main()