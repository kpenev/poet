import numpy as np
import matplotlib.pyplot as plt
import sys
MASS_COL = -2
PERIOD_COL = -4
NUM_BINS = 10

class Star:
    mass = None
    period = None
    def __init__(self, mass, period):
        self.mass = mass
        self.period = period

def mass_compare(x, y):
    if x.mass < y.mass: return -1
    if x.mass > y.mass: return 1
    return 0

if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " filename"
    exit(-1)
filename = sys.argv[1]
masses = []
periods = []
stars = []
with open(filename) as f:
    line = f.readline()
    for line in f:
        elements = line.split()
        mass = float(elements[MASS_COL])
        if mass >= 0.5 and mass <= 1.1:
            stars.append(Star(mass, float(elements[PERIOD_COL])))
            masses.append(mass)
            periods.append(float(elements[PERIOD_COL]))
stars.sort(mass_compare)
bins = np.linspace(0.5, 1.1, NUM_BINS + 1)
print NUM_BINS
for i in range(len(bins)):
    sys.stdout.write(str(bins[i]) + " ")
print ""
curr_bin = 0
for i in range(len(stars)):
    if curr_bin == NUM_BINS - 1 and stars[i].mass >= bins[curr_bin] and \
            stars[i].mass <= bins[curr_bin + 1]:
        #handle special case with last bin (bin edges both inclusive)
        sys.stdout.write(str(stars[i].period) + " ")
        continue
    if stars[i].mass >= bins[curr_bin] and \
            stars[i].mass < bins[curr_bin + 1]:
        sys.stdout.write(str(stars[i].period) + " ")
    else:
        curr_bin += 1
        sys.stdout.write("\n" + str(stars[i].period) + " ")

#for i in range(len(stars)):
 #   print stars[i].mass

fontsize = "x-large"

plt.hist(masses, bins=NUM_BINS)
plt.xlabel("$M_*/M_{sun}$", fontsize=fontsize)
plt.ylabel("# stars", fontsize=fontsize)
plt.figure()

plt.hist(periods, bins=20)
plt.xlabel("Rotation period (days)", fontsize=fontsize)
plt.ylabel("# stars", fontsize=fontsize)
plt.figure()

plt.scatter(masses, periods)
plt.xlabel("$M_*/M_{sun}$", fontsize=fontsize)
plt.ylabel("Rotation period (days)", fontsize=fontsize)
plt.tick_params(labelsize=15)
plt.show()
