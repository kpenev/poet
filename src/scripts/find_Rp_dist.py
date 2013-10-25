import numpy as np
import matplotlib.pyplot as plt
import sys

alpha = -1.92
#number of Earth radii in Jupiter radius
EARTH_TO_JUP = 1.0/11.21
MIN_R = 8*EARTH_TO_JUP
MAX_R = 32*EARTH_TO_JUP
NUM_BINS = 100

def get_prob(R):
    return R**(alpha - 1)

allRadii = np.linspace(MIN_R, MAX_R, NUM_BINS+1)
allProbs = []
midRadii = []
for i in range(NUM_BINS):
    radius = (allRadii[i] + allRadii[i+1])/2.
    midRadii.append(radius)
    allProbs.append(get_prob(radius))
allProbs = np.array(allProbs)
allProbs /= np.sum(allProbs)

print len(allProbs)
for i in range(NUM_BINS + 1):
    sys.stdout.write(str(allRadii[i]) + " ")
print ""
for i in range(NUM_BINS):
    sys.stdout.write(str(allProbs[i]) + " ")
print ""
plt.plot(midRadii, allProbs)
fontsize = "x-large"
plt.xlabel("Planet radius (Jupiter radii)", fontsize=fontsize)
plt.ylabel("Probability density", fontsize=fontsize)
plt.show()
    
