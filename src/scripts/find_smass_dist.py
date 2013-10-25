import csv
import sys
import numpy as np
import matplotlib.pyplot as plt

G = 6.674e-8 #cgs units
solar_radius = 6.955e10 #cm
solar_mass = 1.99e33 #grams
NUM_BINS = 100

if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " filename"
    exit(-1)
filename = sys.argv[1]
masses = []
radii = []
numValid = 0
with open(filename) as csvfile:
    reader = csv.reader(csvfile)
    reader.next()
    reader.next()
    for row in reader:
        if len(row)==0 or row[0]=="" or row[1]=="": continue
        radius = float(row[0])
        log_g = float(row[1])
        radii.append(radius)
        log_mass = log_g - np.log10(G) + 2*np.log10(solar_radius*radius)
        log_mass -= np.log10(solar_mass)
        mass = (10**log_mass)
        if mass >= 0.5 and mass <= 1.1: 
            numValid += 1
            masses.append(mass)

print "Number valid: " + str(numValid)
(data, bins) = np.histogram(masses, bins=NUM_BINS)
data = np.array(data, dtype='float')/np.sum(data)
#print out number of bins and bins themselves
print len(data)
for i in range(len(bins)):
    sys.stdout.write(str(bins[i]) + " ")
print ""
#print out actual data
for i in range(len(data)):
    sys.stdout.write(str(data[i]) + " ")
print ""

plt.hist(masses,bins=NUM_BINS)
#plt.text(fontsize="xx-large")
fontsize = "x-large"
plt.xlabel("$M_*/M_{sun}$", fontsize=fontsize)
plt.ylabel("# stars", fontsize=fontsize)
#plt.title("Truncated mass probability distribution for Kepler targets")
plt.show()
