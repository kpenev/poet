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
#from scipy.integrate import quad
import matplotlib.pyplot as plt
# from matplotlib import gridspec
####

### Declare constants here
#MAX_TERMS = 400
GRID = 10001 # The number of grid points to divide 0<=e<=1 into; should be power of ten plus one

### Classes
#class Classy:

### Define functions
#def function(var1, var2):
#	return var1+var2

#####
# The following are all defining the coefficients
# which we're trying to find expansions of
# I will have removed the poles for all of them eventually some day (1-e**2)**(3/2)*
# For more info, see (documentation url)

# Term inside integral for p0s
def p_0s(u,s,e):
	return (1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))/(1-e*np.cos(u))**2).real

# The next four terms are the constituent parts
# of m=+/-2
def t1(u,s,e):
	return (1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))/(1-e*np.cos(u))**4).real

def t2(u,s,e):
	return (1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.cos(2*u)/(1-e*np.cos(u))**4).real

def t3(u,s,e):
	return (1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.sin(2*u)/(1-e*np.cos(u))**4).real

def t4(u,s,e):
	return (1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.sin(u)/(1-e*np.cos(u))**4).real
	
# The overall coefficient thing
def p_MS(m,s,e):
	
	p0s = integrate.quad(p_0s,0,2*np.pi,args=(s,e))[0]
	
	if abs(m) == 2:
	
		sC1 = 1-e**2    # Some constant
		sC2 = np.sqrt(sC1) # Some other constant
		
		term1 = -sC1*integrate.quad(t1,0,2*np.pi,args=(s,e))[0]
		term2 = sC1*integrate.quad(t2,0,2*np.pi,args=(s,e))[0]
		term3 = 1j*sC2*integrate.quad(t3,0,2*np.pi,args=(s,e))[0]
		term4 = 1j*2*e*sC2*integrate.quad(t4,0,2*np.pi,args=(s,e))[0]
		
		if m > 0:
			term3 = -term3
		elif m < 0:
			term4 = -term4
		
	else:
		term1 = term2 = term3 = term4 = 0
		
	return (1/(2*np.pi))*(p0s+term1+term2+term3+term4).real

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
def main(m=2, s=2, accuracyGoal=3, maxCoeffs=30):
	
	# Sanity check input arguments
		# Correct number of them? Correct type? Acceptable values?
		# If not, complain and exit
		
	# Initialize my variables
	coeffDeg = 0#[0:5] # Current degree(s) of Chebyshev polynomial for which we are finding coefficients
	listOfCoeff = np.zeros((2,1))
	acc = 0
	
	# Optional?
	# Determine range of e we need to deal with (for speed purposes)
		# As s increases, the point where there's noticeable divergence goes further right
		# Obviously need full range for getting coeffs, but maybe not for confirming accuracy?
	
	# Optional.
	# Precalculate p_ms(e) itself to compare to? Can I do that? Would be convenient
	print("hi")
	#Calculate some number of data points for 0<=e<=1
	eList = np.zeros(GRID)
	for i in range(1,GRID):
		eList[i] = i*1/(GRID-1)
	
	# Calculate the value of p_ms for each point on that list
	vec_pms = np.vectorize(p_MS) # Allow the bit I defined to take a vector for arguments, hassle-free
	yList = vec_pms(m,s,eList)
	
	#print(eList)
	#print(yList)
	plt.plot(eList,yList)
	plt.show()
	
	while (acc < accuracyGoal and coeffDeg <= maxCoeffs):
		
		#newTerm = getCoefficients(coeffDeg,eList,yList)
		listOfCoeff = getCoefficients(coeffDeg,eList,yList)
		
		# if does not blow up (use prune, return true or false) # Wait to do this one until you have something running so you can see if it's necessary
		
		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		# double check if you're just doing one coefficient at a time here please
		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		# update: that is not what's happening
		
		#listOfCoeff = np.concatenate((listOfCoeff,np.array([[coeffDeg],[newTerm[0]]])) )
		
		# ListOfCoefficients = prune(listOfC...) # In Mathematica, sometimes we got terms that exploded, but then ones after that
											   #were fine, so I want to be able to remove problematic terms
		
		# accuracy = checkAccuracy(ListOfCoefficients, equation, [range of e], pms) # Check to see how accurate we are
		
		coeffDeg += 1
	
	# Report results
	# Remove the empty first index, unless I change how that works
	print(listOfCoeff)

if __name__ == "__main__":
	main()