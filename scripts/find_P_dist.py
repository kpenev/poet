import numpy as np
import matplotlib.pyplot as plt
import sys

WITH_CUTOFF = True
MIN_PERIOD = 0
MAX_PERIOD = 5.9
NUM_BINS = 100
beta = 0.37
P0 = 1.7
#P0 = 0.68
gamma = 4.1

def get_prob(P):
    if P < MIN_PERIOD or P > MAX_PERIOD: return 0
    if WITH_CUTOFF:
        return (P**(beta - 1))*(1-np.exp(-(P/P0)**gamma))
    return P**(beta - 1)

allP = np.linspace(0, MAX_PERIOD, NUM_BINS + 1)
allProbs = []
midP = []
for i in range(NUM_BINS):
    P = (allP[i] + allP[i+1])/2.0
    midP.append(P)
    allProbs.append(get_prob(P))

allProbs = np.array(allProbs)/np.sum(allProbs)
print NUM_BINS
for i in range(NUM_BINS + 1):
    sys.stdout.write(str(allP[i]) + " ")
print ""
for i in range(NUM_BINS):
    sys.stdout.write(str(allProbs[i]) + " ")
print ""

fontsize = "x-large"
plt.plot(midP, allProbs)
plt.xlabel("Orbital period (days)", fontsize=fontsize)
plt.ylabel("Probability density", fontsize=fontsize)
plt.show()
