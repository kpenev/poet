""" determineHowMuchGridINeedPerUnitGrid.py

Calculates multiple expansions of p_ms coefficients
up to some level of accuracy.
@author: Joshua Aaron Schussler
"""

### Importing relevant libraries
import sys, getopt
import numpy as np
#import scipy.integrate as integrate
import matplotlib.pyplot as plt
#from multiprocessing import Pool
import multiprocessing as mp
import lib_pms as pms
####

def main(m=0, s=1, tolerance=2e-9, kind=3):
    
    print("Starting")
    # Sanity check input arguments
    if (m!=0 and abs(m)!=2 or (not isinstance(m,int))):
        print("Improper m")
        return [[-1]],[-1]
    if (s<0 or (not isinstance(s,int))):
        print("Improper s")
        return [[-1]],[-1]
        
    # Initialize my variables
    maxLoops = 0#30#15#11#np.inf#100
    loops = 0
    goalAchieved = 0
    GRID = 101#11
    zeroEnd = 0
    eZEND = 0
    cutoff=.999
    fNme="calcOutput_multi.txt"
    
    zeeArray = ezList = np.zeros(0)
    
    # Calculate some number of data points for 0<=e<=1
    eList = np.linspace(0,1,GRID)
        
    print("Calculating initial p_ms")
    # Calculate the value of p_ms for each point on that list
    #vec_pms = np.vectorize(p_MS) # Allow the bit I defined to take a vector for arguments, hassle-free
    yList = pms.getCoefficient(m,s,eList,eList.size,kind,fNme,loops-2,tolerance)#vec_pms(m,s,eList)
    print(yList)
    # Identify how long the bit that's just zero is (if it exists)
    print("Getting zero")
    for i in yList:
        if i != 0:
            break
        else:
            zeroEnd = zeroEnd + 1
    print("DEBUG zeroEnd: " + str(zeroEnd))
    # Adjust our setup to account for our new reality (unless nothing changed, in which case that would be a waste of time)
    print("Recalculating")
    if zeroEnd > 1:
        eZEND = eList[zeroEnd-1]
        if eZEND != 1:
            zeeArray = np.zeros(zeroEnd-1)
            ezList = np.linspace(0,1,GRID)[0:zeroEnd-1]
            print("DEBUG Zero is at " + str(eList[zeroEnd]))
            eList = np.linspace(eZEND,1,GRID)
            yList = pms.getCoefficient(m,s,eList,eList.size,kind,fNme,loops-1,tolerance)
        else:
            print("Coefficient is always zero")
            return eZEND,0,cutoff,np.zeros(1)#0
    else:
        print("DEBUG zeroEnd not greater than 1")
    #print(yList)
    print("--- BEGIN EXPLORING ---")
    while ((loops <= maxLoops) and (not goalAchieved)):
        print("Increasing Resolution")
        newGRID = GRID * 2 - 1
        eListTemp = np.linspace(eZEND,1,newGRID)[1::2]
        yListTemp = pms.getCoefficient(m,s,eListTemp,eListTemp.size,kind,fNme,loops,tolerance)
        print("Interpolating")
        interpolant = np.interp(eListTemp,eList,yList)
        
        print("Comparing")
        #goalAchieved = isGoalAchieved(eList,yList,eListTemp,yListTemp)
        print(str(np.max( np.abs(np.array(yListTemp)[eListTemp<=cutoff] - interpolant[eListTemp<=cutoff]) )))
        print(str(eListTemp[np.argmax( np.abs(np.array(yListTemp)[eListTemp<=cutoff] - interpolant[eListTemp<=cutoff]) )]))
        if np.max( np.abs(np.array(yListTemp)[eListTemp<=cutoff] - interpolant[eListTemp<=cutoff]) ) < tolerance:  # do these as vector[0:?]
            goalAchieved = True
        
        if (not goalAchieved):
            print("Combining")
            GRID=newGRID
            newArray = np.zeros(GRID)
            newArray[0::2] = eList
            newArray[1::2] = eListTemp
            eList = newArray
            newArray = np.zeros(GRID)
            newArray[0::2] = yList
            newArray[1::2] = yListTemp
            yList = newArray
        
        loops += 1
    
    print("--- --------------- ---")
    print("Finishing")
    # Finish
    #feL = np.append(ezList,eList)
    #fyL = np.append(zeeArray,yList)
    # Report results
    #print(ezList)
    #print(eList)
    print("We used this many grid points: " + str(GRID))
    #plt.plot(feL,fyL,'o')
    #plt.show()
    #plt.plot(ezList,zeeArray,'o')
    #plt.plot(feL,fyL,'+')
    #plt.show()
    #plt.semilogy(eList,yList,'o')
    #plt.semilogy(eList,-1*yList,'+')
    #plt.show()
    print(goalAchieved)
    
    return eZEND,GRID,cutoff,yList

if __name__ == "__main__":
    mp.set_start_method('spawn')
    main()