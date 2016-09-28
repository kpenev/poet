import numpy as np
import matplotlib.pyplot as plt
import sys

AGE_COL = 0
SPIN_COL = 1
SEMI_COL = 2
NUM_ELEMENTS=5

class Point:
    age = None
    semi = None
    spin = None
    def __init__(self, age, semi, spin):
        self.age = age
        self.semi = semi
        self.spin = spin

def age_compare(x, y):
    if x.age < y.age: return -1
    if x.age > y.age: return 1
    return 0

if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " filename"
    exit(-1)
filename = sys.argv[1]
f = open(filename)
ages = []
spins = []
semis = []
line = f.readline() #get rid of first 2 lines
line = f.readline()
line = f.readline()
datapoints = []
while line != "":
    elements = line.split()
    if len(elements) != NUM_ELEMENTS: break
    age = float(elements[AGE_COL])
    semi = float(elements[SEMI_COL])
    spin = float(elements[SPIN_COL])
    datapoints.append(Point(age, semi, spin))
    line = f.readline()
datapoints.sort(age_compare)
ages = []
spins = []
semis = []
for i in range(len(datapoints)):
    ages.append(datapoints[i].age)
    spins.append(datapoints[i].spin)
    if datapoints[i].semi == 1:
        semis.append(np.nan)
    else: semis.append(100*datapoints[i].semi)
plt.plot(ages, semis, label="a (0.01 AU)")
plt.xlabel("Age (Gyr)")
plt.plot(ages, spins, label="$\omega_{conv}$ (rad/day)")
plt.legend(frameon=False)
plt.show()
