""" Some initial documentation

(details)
@author: Joshua Aaron Schussler
"""

### Import things here
# import os
# import math
# import numpy as np
# import scipy.stats as stats
import scipy.integrate as integrate
# import matplotlib.pyplot as plt
# from matplotlib import gridspec
####

### Declare constants here
#MAX_TERMS = 400

### Classes
#class Classy:

### Define functions
#def function(var1, var2):
#	return var1+var2

# The overall coefficient thing
def pMS(m,s):
	
	h

# Documentation goes here
# Gets the specified coefficient of the (Chebyshev?) expansion
def getCoefficients(equation, [termsToFind], ):
	for each term to find
		result[] += findChebyshevExpansion(equation, termsToFind[i]) # There will be some built in thing, I'm sure
																	 # If not, I guess I'll be defining this function later
	return result

# Removes junk terms
def prune(lOC):
	for each term in lOC
		if term is greater than 10
			delete the term from lOC
	return lOC

# Checks if we're accurate enough
def checkAccuracy(ListOfCoefficients, equation, [range of e], pms):
	# pms is a numerically integrated version of p_ms that we did beforehand for just this comparison
	I'm not entirely sure how to do this
		1. Choose regularly spaced positions in range and get the percent error?
		2. is there some other way to quantify how closely one function emulates another?
	if we have multiple terms, do we average them for the overall accuracy or combine them in some other way?
	else if there's only one term, problem solved I guess
	return accuracy

### Run stuff
#function(0,1)

--start execution
# I'm not sure if I want the user to specify a range of m and s or just a specific p_ms (the latter would need
# the program to be run individually for each p_ms coefficient)
def main function(m=0, s=1, accuracyGoal, maxCoeffs):
	
	Sanity check input arguments
		Correct number of them? Correct type? Acceptable values?
		If not, complain and exit
	
	Set up the exact equation we're going to use
		if m = 0, all terms but first shouldn't be there
		if m=+/- 2, those terms should be there
		Options:
			1. Variable referring to equation
				a. Set some constant that makes the equation be the form I want (possibly slow, might spend time recalculating an integral
				 ... only to multiply by zero)
				b. Change which equation the variable points to (pointers!)
			2. ??? how does numpy work w/ numerical integrating equations?
		
	Determine range of e we need to deal with (for speed purposes)
		As s increases, the point where there's noticeable divergence goes further right
		Obviously need full range for getting coeffs, but maybe not for confirming accuracy?
		
	Precalculate p_ms(e) itself to compare to? Can I do that? Would be convenient
	
	loop while not accurate enough or number of terms is below maxCoeffs (?)
		
		ListOfCoefficients += [ [termsToFind] , getCoefficients(equation, [termsToFind], ) ] # Find however many more coefficients
																							 # I think I want lOC to be [ChebNumber,Coef] (2d)
																							 # That way the next step doesn't cause problems
		
		ListOfCoefficients = prune(listOfC...) # In Mathematica, sometimes we got terms that exploded, but then ones after that
											   # were fine, so I want to be able to remove problematic terms
		
		accuracy = checkAccuracy(ListOfCoefficients, equation, [range of e], pms) # Check to see how accurate we are
	
	Report results
