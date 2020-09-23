""" Some initial documentation

(details)
@author: Joshua Aaron Schussler
"""

### Import things here
import sys, getopt
# import os
# import math
import numpy as np
# import scipy.stats as stats
import scipy.integrate as integrate
from integrate import quad
# import matplotlib.pyplot as plt
# from matplotlib import gridspec
####

### Declare constants here
#MAX_TERMS = 400
GRID = 101 # The number of grid points to divide 0<=e<=1 into; should be power of ten plus one

### Classes
#class Classy:

### Define functions
#def function(var1, var2):
#	return var1+var2

#####
# The following are all defining the coefficients
# which we're trying to find expansions of
# For more info, see (documentation url)

# Term inside integral for p0s
def p_0s(x,m,s,e):
	return exp(i*s*(x-e*sin(x)))/(1-e*cos(u))**2

# The next four terms are the constituent parts
# of m=+/-2
def t1(x,m,s,e):
	return exp(i*s*(x-e*sin(x)))/(1-e*cos(u))**4

def t2(x,m,s,e):
	return exp(i*s*(x-e*sin(x)))*cos(2*x)/(1-e*cos(u))**4

def t3(x,m,s,e):
	return exp(i*s*(x-e*sin(x)))*sin(2*x)/(1-e*cos(u))**4

def t4(x,m,s,e):
	return exp(i*s*(x-e*sin(x)))*sin(x)/(1-e*cos(u))**4
	
# The overall coefficient thing
def p_MS(m,s,e):
	
	p0s = quad(p_0s,0,2*pi,args=(m,s,e))
	
	if abs(m) == 2:
	
		sC1 = 1-e**2    # Some constant
		sC2 = sqrt(sC1) # Some other constant
		
		term1 = -sC1*quad(t1,0,2*pi,args=(m,s,e))
		term2 = sC1*quad(t2,0,2*pi,args=(m,s,e))
		term3 = i*sC2*quad(t3,0,2*pi,args=(m,s,e))
		term4 = i*2*e*sC2*quad(t4,0,2*pi,args=(m,s,e))
		
		if m > 0:
			term3 = -term3
		elif m < 0:
			term4 = -term4
		
	else:
		term1=term2=term3=term4 = 0
		
	return (1/(2*pi))*(p0s+term1+term2+term3+term4)

# Documentation goes here
# Gets the specified coefficient of the (Chebyshev?) expansion
def getCoefficients(coeffDeg,eList,yList):
	#for each term to find
		#result[] += findChebyshevExpansion(equation, termsToFind[i]) # There will be some built in thing, I'm sure
																	 #If not, I guess I'll be defining this function later
	
	return np.polynomial.chebyshev.chebfit(eList,yList,coeffDeg)

# Removes junk terms
# def prune(lOC):
	# for each term in lOC
		# if term is greater than 10
			# delete the term from lOC
	# return lOC

# Checks if we're accurate enough
#def checkAccuracy(ListOfCoefficients, equation, [range of e], pms):
	# pms is a numerically integrated version of p_ms that we did beforehand for just this comparison
#	I'm not entirely sure how to do this
#		1. Choose regularly spaced positions in range and get the percent error?
#		2. is there some other way to quantify how closely one function emulates another?
#	if we have multiple terms, do we average them for the overall accuracy or combine them in some other way?
#	else if there's only one term, problem solved I guess
#	return accuracy

# I'm not sure if I want the user to specify a range of m and s or just a specific p_ms (the latter would need
# the program to be run individually for each p_ms coefficient)
def main(m=0, s=1, accuracyGoal, maxCoeffs):
	
	# Sanity check input arguments
		# Correct number of them? Correct type? Acceptable values?
		# If not, complain and exit
		
	# Initialize my variables
	coeffDeg = [0:5] # Current degree(s) of Chebyshev polynomial for which we are finding coefficients
	listOfCoeff = np.zeros(2,however many coefficients)
	
	# Optional?
	# Determine range of e we need to deal with (for speed purposes)
		# As s increases, the point where there's noticeable divergence goes further right
		# Obviously need full range for getting coeffs, but maybe not for confirming accuracy?
	
	# Optional.
	# Precalculate p_ms(e) itself to compare to? Can I do that? Would be convenient
	
	#Calculate some number of data points for 0<=e<=1
	eList = np.zeros(GRID)
	for i in range(1,GRID-1):
		eList[i] = i*1/(GRID-1)
	
	# Calculate the value of p_ms for each point on that list
	vec_pms = vecorize(p_MS) # Allow the bit I defined to take a vector for arguments, hassle-free
	yList = vec_pms(m,s,eList)
	
	while (acc < accuracyGoal && length(listOfCoeff[1]) < maxCoeffs):
		
		newTerms = getCoefficients(termsToFind,eList,yList)
		listOfCoeff = np.concatenate((listOfCoeff,np.array([termsToFind,newTerms])) )
		
		# Wait to do this one until you have something running so you can see if it's necessary
		# ListOfCoefficients = prune(listOfC...) # In Mathematica, sometimes we got terms that exploded, but then ones after that
											   #were fine, so I want to be able to remove problematic terms
		
		# accuracy = checkAccuracy(ListOfCoefficients, equation, [range of e], pms) # Check to see how accurate we are
	
	# Report results

if __name__ == "__main__":
	main()