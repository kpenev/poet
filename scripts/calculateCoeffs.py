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
	
	##### ~~!!!!!!!!~
	### How good is this calculation??? What can I do about improving that as needed?
	
	limitBreak = 10*s
	if limitBreak < 50:
		limitBreak = 50
	
	if e == 1: # Special case of we-know-what-this-should-be
	
		if m == 0:
			return 1
		else:
			return 0
	
	else:
		p0s = integrate.quad(p_0s,0,2*np.pi,args=(s,e),limit=limitBreak)[0] ######### INCREASE LIMIT
															#### literally directly from s
															### limit: at least 10*s
															# don't go below default 50
		
		if abs(m) == 2:
		
			sC1 = 1-e**2    # Some constant
			sC2 = np.sqrt(sC1) # Some other constant
			
			term1 = -sC1*integrate.quad(t1,0,2*np.pi,args=(s,e),limit=limitBreak)[0]
			term2 = sC1*integrate.quad(t2,0,2*np.pi,args=(s,e),limit=limitBreak)[0]
			term3 = 1j*sC2*integrate.quad(t3,0,2*np.pi,args=(s,e),limit=limitBreak)[0]
			term4 = 1j*2*e*sC2*integrate.quad(t4,0,2*np.pi,args=(s,e),limit=limitBreak)[0]
			
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
	
	#print(np.polynomial.chebyshev.chebfit(eList,yList,coeffDeg,True)[0])
	#print(np.polynomial.chebyshev.chebfit(eList,yList,coeffDeg,True)[1])
	#bob, george = np.polynomial.chebyshev.chebfit(eList,yList,coeffDeg,None,True)
	#print(eList.size)
	#print(eList.ndim)
	#print((eList.T).size)
	#print((eList.T).ndim)
	
	# 2*e-1
	
	return np.polynomial.chebyshev.chebfit(2*eList.T-1,yList.T,coeffDeg,None,True)

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
def main(m=0, s=10, accuracyGoal=3, maxCoeffs=30):
	
	# Sanity check input arguments
		# Correct number of them? Correct type? Acceptable values?
		# If not, complain and exit
		
	# Initialize my variables
	# Variables related to looping
	coeffDeg = 0# Current degree(s) of Chebyshev polynomial for which we are finding coefficients
	notDone = 1 # Whether we're done looping because we've found a good result or not
	# Variables related to accuracy
	acc = 10000000
	#bestResult = fsf # residual,{other relevant bits}
	resid = np.array(100) # Residual
	# Variables related to reporting
	listOfCoeff = np.zeros((2,1))
	# Miscellaneous variables
	startAt = 0 # A number between 0 and 1 at which interesting things start happening
	
	# Determine range of e we need to deal with (for speed purposes)
		# As s increases, the point where there's noticeable divergence goes further right
		# Obviously need full range for getting coeffs, but maybe not for confirming accuracy?
		
		# Do guess your number thing? bouncing back and forth and stuff
	
	# Calculate some number of data points for 0<=e<=1
	# Make two; the first is for getting coefficients, the second is for checking the accuracy of them
	
	##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	## LINSPACE <<<<<< --------- !!!!!!
	##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	# ^
	# |
	# |
	# | (but need transformations for big s)
	# |
	# |
	
	if 0:
		eList = np.zeros(GRID)
		for i in range(np.around(startAt*GRID),GRID):
			eList[i] = i*(1-startAt)/(GRID-1)
	eList = np.linspace(startAt,1,GRID)
	#eListCheck = np.zeros(GRID*2-1)
	#for i in range(np.around(startAt*(GRID*2-1)),GRID*2-1):
	#	eListCheck[i] = i*(1-startAt)/(GRID*2-2)
		
	#####~~!!!! ^^^ - dynamic resolution for when there's only a hundred points doing interesting things? evidence
	### calculating a middle point and then comparing to fit
	### or just two res and res/2 and see if residual changes much
	
	# Calculate the value of p_ms for each point on that list
	vec_pms = np.vectorize(p_MS) # Allow the bit I defined to take a vector for arguments, hassle-free
	yList = vec_pms(m,s,eList)
	
	#print(eList)
	#print(yList)
	
	# This bit is plotting
	#plt.plot(eList,yList)
	#plt.show()
	
	# while overall
		# loop number coefficients
			# loop residual
	# Find point where residual was lowest?
	# Re-get that point
	# Wait why do I need to do that? How many things do I really have available to adjust?
	
	#while notDone:
	
	#while coeffDeg <= maxCoeffs:
	
	while (resid.size and coeffDeg <= maxCoeffs): ####### <<---- residual (maybe make a plot of precision vs. (degree/resolution)
		
		#newTerm = getCoefficients(coeffDeg,eList,yList)
		listOfCoeff, diag = getCoefficients(coeffDeg,eList,yList) #########<---------- I may need to flip the points to go -1 to 1
		print(diag[0])
		print(diag[0].size)
		resid = diag[0]
		# if does not blow up (use prune, return true or false) # Wait to do this one until you have something running so you can see if it's necessary
		
		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		# double check if you're just doing one coefficient at a time here please
		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		# update: that is not what's happening
		
		#listOfCoeff = np.concatenate((listOfCoeff,np.array([[coeffDeg],[newTerm[0]]])) )
		
		# accuracy = checkAccuracy(ListOfCoefficients, equation, [range of e], pms) # Check to see how accurate we are
		####### accuracy: residuals < some cutoff value OR all results for second array less than some percent error?
		if abs((acc - diag[0])) < 1:
			
			# We're accurate enough
			maxCoeffs = coeffDeg
		
		acc = diag[0]
		coeffDeg += 1
	
	# The loop is going to increase coeffDeg too far, so fix that
	coeffDeg = coeffDeg - 1
	#f
	
	# Report results
	# Remove the empty first index, unless I change how that works
	#print(listOfCoeff)
	#print(diag)
	#print(diag[0])
	#print(diag[1])
	#print(diag[2])
	
	theCheb = np.polynomial.chebyshev.Chebyshev(listOfCoeff)
	plt.plot(eList,yList,'o')
	#dirtX = theCheb.linspace(GRID)[0]
	#dirtX = dirtX[int(np.around(dirtX.size/2)):int(dirtX.size)]
	#dirtX = np.take(dirtX,eList[(GRID-1)/2:GRID]*(GRID-1))
	dirtY = np.polynomial.chebyshev.chebval(2*eList-1,listOfCoeff)
	#dirtY = theCheb.linspace(GRID)[1]
	#dirtY = dirtY[int(np.around(dirtY.size/2)):int(dirtY.size)]
	plt.plot(eList,dirtY)
	plt.show()
	
	#sql alchemy

if __name__ == "__main__":
	main()