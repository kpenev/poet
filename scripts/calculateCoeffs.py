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

# The next four terms are the constituent parts of m=+/-2
#  t3 and t4 are later multiplied by 1j; because quad can't handle
#  complex numbers and we ultimately only want the real part, we
#  handle the flipping and subsequent negative here
def t1(u,s,e):
	return ((1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))/(1-e*np.cos(u))**4)).real

def t2(u,s,e):
	return (1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.cos(2*u)/(1-e*np.cos(u))**4).real

def t3(u,s,e):
	return -((1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.sin(2*u)/(1-e*np.cos(u))**4)).imag

def t4(u,s,e):
	return -((1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.sin(u)/(1-e*np.cos(u))**4)).imag
	
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
			term3 = sC2*integrate.quad(t3,0,2*np.pi,args=(s,e),limit=limitBreak)[0]
			term4 = 2*e*sC2*integrate.quad(t4,0,2*np.pi,args=(s,e),limit=limitBreak)[0]
			
			if m > 0:
				term3 = -term3
			elif m < 0:
				term4 = -term4
			
		else:
			term1 = term2 = term3 = term4 = 0
		
	return (1/(2*np.pi))*(p0s+term1+term2+term3+term4)

# Gets the specified coefficient of the (Chebyshev?) expansion
def getCoefficients(coeffDeg,eList,yList):
	
	return np.polynomial.chebyshev.chebfit(2*eList.T-1,yList.T,coeffDeg,None,True)

def main(m=2, s=10, accuracyGoal=1, maxCoeffs=np.inf):#1e-6
	
	print("Starting")
	# Sanity check input arguments
	if (m!=0 and abs(m)!=2 or (not isinstance(m,int))):
		print("Improper m")
		return [[-1]],[-1]
	if (s<0 or (not isinstance(s,int))):
		print("Improper s")
		return [[-1]],[-1]
	#more
		
	# Initialize my variables
	coeffDeg = 0 # Current degree(s) of Chebyshev polynomial for which we are finding coefficients
	resid = [] # Residual
	listOfCoeff=[]
	
	# Calculate some number of data points for 0<=e<=1
	eList = np.linspace(0,1,GRID)
		
	print("Calculating p_ms")
	# Calculate the value of p_ms for each point on that list
	vec_pms = np.vectorize(p_MS) # Allow the bit I defined to take a vector for arguments, hassle-free
	yList = vec_pms(m,s,eList)
	
	print("Fitting")
	while coeffDeg <= maxCoeffs:
		
		newCoeff, diag = getCoefficients(coeffDeg,eList,yList)
		listOfCoeff.append(newCoeff)
		resid.append(diag[0][0]) # Diag is an array of four items, the one we want is itself an array of a single element
		
		if resid[coeffDeg] < accuracyGoal: # Don't change this below 1e-7
			
			# We're accurate enough
			maxCoeffs = coeffDeg
		
		coeffDeg += 1
	
	# The loop is going to increase coeffDeg too far, so fix that
	coeffDeg = coeffDeg - 1
	
	# Report results
	print("We used this many coefficients: " + str(coeffDeg))
	theCheb = np.polynomial.chebyshev.Chebyshev(listOfCoeff[coeffDeg])
	plt.plot(eList,yList,'o')
	dirtY = np.polynomial.chebyshev.chebval(2*eList-1,listOfCoeff[coeffDeg])
	plt.plot(eList,dirtY)
	plt.show()
	
	return listOfCoeff, resid

if __name__ == "__main__":
	main()