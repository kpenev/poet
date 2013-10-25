import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
NUM_BINS = 100
MASS_COL = 0
RADIUS_COL = 1
if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " filename"
    exit(-1)
filename = sys.argv[1]
masses = []
radii = []
with open(filename) as csvfile:
    reader = csv.reader(csvfile)
    reader.next()
    reader.next()
    for row in reader:
        if len(row)==0 or row[0]=="" or row[1]=="": continue
        mass = float(row[MASS_COL])
        radius = float(row[RADIUS_COL])
        masses.append(mass)
        radii.append(radius)

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

fontsize = "x-large"

plt.hist(masses, bins=NUM_BINS)
plt.xlabel("Planet mass (Jupiter masses)", fontsize=fontsize)
plt.ylabel("# planets", fontsize=fontsize)
plt.figure()

plt.scatter(radii, masses)
plt.ylabel("Planet mass (Jupiter masses)", fontsize=fontsize)
plt.xlabel("Planet radius (Jupiter radii)", fontsize=fontsize)
plt.show()
