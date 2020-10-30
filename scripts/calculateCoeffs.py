""" calculateCoeffs.py

Calculates multiple expansions of p_ms coefficients
up to some level of accuracy.
@author: Joshua Aaron Schussler
"""

### Importing relevant libraries
import sys, getopt
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
####

### Constants
GRID = 10001 # The number of grid points to divide 0<=e<=1 into; should be power of ten plus one

#####
# The following are all defining the coefficients
# which we're trying to find expansions of
# I will have removed the poles for all of them eventually some day (1-e**2)**(3/2)*
# For more info, see https://kpenev.github.io/poet/inclination_eccentricity_pms1.html

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
	
# The p_ms coefficient proper
def p_MS(m,s,e):
	
	# We want to adjust the numerical integration limit as s
	# increases, but we don't want to go below the default value
	# for small s
	limitBreak = 10*s
	if limitBreak < 50:
		limitBreak = 50
	
	if e == 1: # Special case of we-know-what-this-should-be
	
		if m == 0:
			return 1
		else:
			return 0
	
	else:
		p0s = integrate.quad(p_0s,0,2*np.pi,args=(s,e),limit=limitBreak)[0]
		
		if abs(m) == 2:
		
			# Define a couple of constants that show up across terms
			sC1 = 1-e**2
			sC2 = np.sqrt(sC1)
			
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

# Gets the specified coefficient of the (Chebyshev?) expansion
def getCoefficients(coeffDeg,eList,yList):
	
	return np.polynomial.chebyshev.chebfit(2*eList.T-1,yList.T,coeffDeg,None,True)

# Checks if we're accurate enough
#def checkAccuracy(ListOfCoefficients, equation, [range of e], pms):
#	return accuracy

# I'm not sure if I want the user to specify a range of m and s or just a specific p_ms (the latter would need
# the program to be run individually for each p_ms coefficient)
def main(m=2, s=100, accuracyGoal=3, maxCoeffs=np.inf):
	print("Starting\n")
	# Sanity check input arguments
		# Correct number of them? Correct type? Acceptable values?
		# If not, complain and exit
		
	# Initialize my variables
	# Variables related to looping
	coeffDeg = 0# Current degree(s) of Chebyshev polynomial for which we are finding coefficients
	notDone = 1 # Whether we're done looping because we've found a good result or not
	# Variables related to accuracy
	acc = 10000000
	resid = np.array(2) # Residual
	# Variables related to reporting
	#listOfCoeff = np.zeros((2,1))
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
	#eList = np.linspace(.95,1,501)
	#eList = np.array([0.999,0.9991,0.9992,0.9993,0.9994,0.9995,0.9996,0.9997,0.9998,0.9999,1])
	#eListCheck = np.zeros(GRID*2-1)
	#for i in range(np.around(startAt*(GRID*2-1)),GRID*2-1):
	#	eListCheck[i] = i*(1-startAt)/(GRID*2-2)
		
	#####~~!!!! ^^^ - dynamic resolution for when there's only a hundred points doing interesting things? evidence
	### calculating a middle point and then comparing to fit
	### or just two res and res/2 and see if residual changes much
	print("Calculating p_ms\n")
	# Calculate the value of p_ms for each point on that list
	vec_pms = np.vectorize(p_MS) # Allow the bit I defined to take a vector for arguments, hassle-free
	yList = vec_pms(m,s,eList)
	
	# while overall
		# loop number coefficients
			# loop residual
	# Find point where residual was lowest?
	# Re-get that point
	# Wait why do I need to do that? How many things do I really have available to adjust?
	
	#while notDone:
	
	#while coeffDeg <= maxCoeffs:
	print("Fitting\n")
	while (resid.size and coeffDeg <= maxCoeffs): ####### <<---- residual (maybe make a plot of precision vs. (degree/resolution)
		
		#newTerm = getCoefficients(coeffDeg,eList,yList)
		listOfCoeff, diag = getCoefficients(coeffDeg,eList,yList) #########<---------- I may need to flip the points to go -1 to 1
		print(diag[0])
		print(diag[0].size)
		resid = diag[0]
		
		if diag[0] < 1e-6: # Don't change this below 1e-7
			
			# We're accurate enough
			maxCoeffs = coeffDeg
		
		acc = diag[0]
		coeffDeg += 1
	
	# The loop is going to increase coeffDeg too far, so fix that
	coeffDeg = coeffDeg - 1
	
	# Report results
	try:
		print("We used this many coefficients: " + str(coeffDef))
	except:
		print("We used this many coefficients: ")
		print(coeffDeg)
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
	## database normal forms
	# relational databases
	# read about sql

if __name__ == "__main__":
	main()