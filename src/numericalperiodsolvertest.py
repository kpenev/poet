import numpy as np
import scipy
import subprocess
import time
import sys
import os
import signal



#can't exactly be pi or the simulation cancels out
incl = np.linspace(0, np.pi * 0.98, 10)
logQinr = np.array([4, 5, 6])
logQ = np.array([6, 7, 8])

TIME_TO_WAIT = 120

planetdata = np.loadtxt('initialplanetsearch4.txt', skiprows=4, delimiter=',', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12), \
             dtype={'names':('NAMES', 'PMASS', 'STARMASS', 'PRADIUS', 'A', 'ORBPERIOD', 'INCL', 'OBLI', 'OBLIMINUS', 'OBLIPLUS', 'AGE', 'AGEMINUS', 'AGEPLUS'), \
             'formats':('S10','f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')})

for i in range(planetdata.shape[0]):
    planetdata['NAMES'][i] = planetdata['NAMES'][i].rstrip()
    planetdata['NAMES'][i] = planetdata['NAMES'][i].replace(' ', '')


def readEvolution(location):
    #output columns will be t, period, a, convincl, radincl, Lconv, Lrad, Iconv, Irad, I, mode
    data = np.genfromtxt(location, skiprows=1, usecols=(0, 1, 2, 3, 4, 10), dtype={'names': ('t', 'ORBPERIOD', 'a', 'convincl', 'radincl', 'mode'), 'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'S15')})
    return data


#array should contain the age, period, Pmass, StarMass, initial obliquity, and the rotational parameters
#returns tuple with string for the command and the string for the output file
#the outf only depends on the input, not the output period and its iteration
#solverEnd is either '+' or '-' or 'c' to indicate how to save the file as bk iteratio or ak iteratio or c interation
def buildCommandString(array, period, solverEnd):
    solverEnd = solverEnd
    initialentry = './poet' + ' '
    inputinfo = '--Mstar=' + str(array['STARMASS'][0]) + ' ' + '--Mplanet=' + str(array['PMASS'][0]) + ' ' \
            + '--init-orbit-period=' + str(period) + ' ' + '--require-ages=' + str(array['AGE'][0]) + ' '\
            + '--init-inclination=' + str(array['OBLI'][0]) + ' ' + '--lgQ=' + str(array['lgQ'][0]) + ' '\
            + '--lgQinr=' + str(array['lgQinr'][0]) + ' ' + '--p-disk=' + str(array['p-disk'][0]) + ' '\
            + '-K ' + str(array['K'][0]) + ' ' + '--core-env-coupling-timescale=' + str(array['tcoup'][0]) + ' '\
            + '--high-mass-wind-sat-w=' + str(array['wsat'][0]) + ' '+ '--t-disk=' + str(array['tdisk'][0]) + ' ' + '--tmax=' + str(1.01 * array['AGE'][0]) + ' '
    
    outputinfo = '--output-columns=t,Porb,a,convincl,radincl,Lconv,Lrad,Iconv,Irad,I,mode '
    location = './solver' + solverEnd + '.evol'
    locationinfo = '--output=' + location
    process = subprocess.Popen(str.split(initialentry + inputinfo + outputinfo + locationinfo, ' '))
    
    t_zero = time.time()
    seconds_passed = 0
    
    while (process.poll() is None and seconds_passed < TIME_TO_WAIT):
        seconds_passed = time.time()-t_zero

    if seconds_passed > TIME_TO_WAIT:
        print 'Time out'
        os.kill(process.pid, signal.SIGKILL)
        return (str.split(initialentry + inputinfo + outputinfo + locationinfo, ' '), -1)
        
    
        
    return (str.split(initialentry + inputinfo + outputinfo + locationinfo, ' '), location)
    


#given a planet
#no need to instantiate all the simulations
#this uses the algorithm of Dekker's method
#we interpolate around the currently measured period instead of zero
#returns an array with initial period and obliquity at that age interpolated to
#both are negative -1 if no solution is found
def solve(array):

    solution = np.array([(0, 0)], dtype=np.dtype([('initialperiod', 'f8'), ('OBLI','f8')]))
    
    UPPER = '+'
    LOWER = '-'

    age = array['AGE'][0]
    #the following are $b_1$ and $a_1$ as in the notation of Dekker's method
    periodUpper = 10 
    
    
    periodLower = 0.2 * array['ORBPERIOD'][0]
    period = array['ORBPERIOD'][0]

    upperProcessString, upperOutf = buildCommandString(array, periodUpper, UPPER)
    

    lowerProcessString, lowerOutf = buildCommandString(array, periodLower, LOWER)
 
    
    #wait for simulations to finish before continuing to work with them
    #exit = [p.wait() for p in process1, process2]

    datatempUpper = readEvolution(upperOutf)
    datatempLower = readEvolution(lowerOutf)

    maskindicesupper = np.isclose(datatempUpper['t'], age)
    maskindiceslower = np.isclose(datatempLower['t'], age)


    #didn't get to age
    #these simulations were never finishing
    if (not np.any(maskindicesupper)):
        while not np.any(maskindicesupper):
            solution['initialperiod'] = -1
            solution['OBLI'] = -1
            return solution
                
    #the estimates are the output f(b) and f(a) respectively, but notice that they are offset from zero
    #the array call is needed as sometimes it returns two
    periodUpperEstimate = np.mean(datatempUpper['ORBPERIOD'][maskindicesupper])
    print periodUpperEstimate
    
    if (not np.any(maskindiceslower)):
        periodLowerEstimate = 0
    else:
        periodLowerEstimate = np.mean(datatempLower['ORBPERIOD'][maskindiceslower])
    #make sure that solution is inside region
    
    if np.isnan(periodLowerEstimate):
        periodLowerEstimate = 0
    
    if np.isnan(periodUpperEstimate):
        count = 0
        while np.isnan(periodUpperEstimate) and count < 2:
            count = count + 1
            print 'Loop 5'
            periodUpper = periodUpper * 1.4
            upperProcessString, upperOutf = buildCommandString(array, periodUpper, UPPER)
            
            if upperOutf == -1:
                solution['initialperiod'] = -1
                solution['OBLI'] = -1
                return solution 
            datatempUpper = readEvolution(upperOutf)
            maskindicesupper = np.isclose(datatempUpper['t'], age)
            
            if (not np.any(maskindicesupper)):
                solution['initialperiod'] = -1
                solution['OBLI'] = -1
                return solution
            
            periodUpperEstimate = np.mean(datatempUpper['ORBPERIOD'][maskindicesupper])
    
    #couldn't bound top
    
    if np.isnan(periodUpperEstimate):
        periodUpperEstimate=0
        
    while (periodUpperEstimate - period) < 0:
        periodUpper = periodUpper * 1.4
        upperProcessString, upperOutf = buildCommandString(array, periodUpper, UPPER)
        
        if upperOutf == -1:
            print 'Solution 1'
            solution['initialperiod'] = -1
            solution['OBLI'] = -1
            return solution
        datatempUpper = readEvolution(upperOutf)
        maskindicesupper = np.isclose(datatempUpper['t'],age)
        #solution wasn't bounded for Hat p 17 b
        #when the step size was increased
        if (not np.any(maskindicesupper)):
            print 'SOlution 2'
            solution['initialperiod'] = -1
            solution['OBLI'] = -1
            return solution
        periodUpperEstimate = np.mean(datatempUpper['ORBPERIOD'][maskindicesupper])

    while (periodLowerEstimate - period) > 0:
        periodLower = periodLower * 0.5
        lowerProcessString, lowerOutf = buildCommandString(array, periodLower, LOWER)
        
        if lowerOutf == -1:
            periodLowerEstimate = 0
        else:
            datatempLower = readEvolution(lowerOutf)
            maskindiceslower = np.isclose(datatempLower['t'],age)
            #solution wasn't bounded for Hat p 17 b
            #when the step size was increased
            if (not np.any(maskindiceslower)):
                periodLowerEstimate = 0
            periodLowerEstimate = np.mean(datatempLower['ORBPERIOD'][maskindiceslower])


    fractionalError = 1.0 * np.abs(period-periodUpperEstimate)/period
    #just larger to better iterate, we want .1% difference so at least iterate if the %error changes
    
    
    bk = periodUpper
    ak = periodLower
    
    if ((periodUpperEstimate-period) * (periodLowerEstimate - period)) > 0:
            solution['initialperiod'] = -1
            solution['OBLI'] = -1
            return solution

    if (np.abs(periodLowerEstimate-period) < np.abs(periodUpperEstimate - period)):
        temp = periodLowerEstimate
        periodLowerEstimate = periodUpperEstimate
        periodUpperEstimate = temp

        temp2 = ak
        ak = bk
        bk = temp2

    c = ak
    periodEstimateC = periodLowerEstimate
    
    #use of d in the solver
    mflag = True 
    
    #set d to get rid of python error, but d won't be used until after the first iteration because of the mflag
    #convergence of the guess points to interpolate between
    tolerancedelta = 0.01
    fractionalErrorChange = 1
    while fractionalError > 0.001 and (np.abs(ak - bk) > 0.001):
        #inverse quadratic interpolation step
        s = 0
        if (periodEstimateC != periodLowerEstimate) and (periodUpperEstimate != periodEstimateC):
            s = ak * (periodUpperEstimate - period) * (periodEstimateC - period) / ((periodLowerEstimate - periodUpperEstimate) * (periodLowerEstimate - periodEstimateC)) + \
                bk * (periodLowerEstimate - period) * (periodEstimateC - period) / ((periodUpperEstimate - periodLowerEstimate) * (periodUpperEstimate - periodEstimateC)) + \
                c * (periodLowerEstimate - period) * (periodUpperEstimate - period) / ((periodEstimateC - periodLowerEstimate) * (periodEstimateC - periodUpperEstimate))


        else:
            s = bk - 1.0 * (periodUpperEstimate - period) * (bk - ak) / (periodUpperEstimate - periodLowerEstimate)

        if (not (s >= ((3.0 * ak + bk)/4.0) and (s <= bk)) or (mflag and (np.abs(s-bk) >= (np.abs(bk-c)/2.0))) or \
                (not mflag and (np.abs(s-bk) >= (np.abs(c-d)/2.0))) or (mflag and np.abs(bk -c) < np.abs(tolerancedelta)) or \
                (not mflag and np.abs(c - d) < np.abs(tolerancedelta))):
            s = (ak + bk)/2.0
            mflag = True
        else:
            mflag = False

        #need to calculate the period that s iterates to
        sProcessString, sOutf = buildCommandString(array, s, 's')
        if sOutf == -1:
            solution['initialperiod'] = -1
            solution['OBLI'] = -1
            return solution
        
        datatempS = readEvolution(sOutf)
        maskindicesS = np.isclose(datatempS['t'], age)
        if (not np.any(maskindicesS)):
            periodEstimateS = 0
        else:
            periodEstimateS = np.mean(datatempS['ORBPERIOD'][maskindicesS])
        
        #should never happen
        if (np.isnan(periodEstimateS)):
            periodEstimateS = 0
        
        
        d = c
        c = bk
        periodEstimateC = periodUpperEstimate
        

        if ((periodLowerEstimate - period) * (periodEstimateS - period)) < 0:
            print 'Changing b'
            bk = s
            periodUpperEstimate = periodEstimateS
        else:
            print 'Changing a'
            ak = s
            periodLowerEstimate = periodEstimateS

        if (np.abs(periodLowerEstimate-period) < np.abs(periodUpperEstimate - period)):
            print 'Switching a and b'
            temp = periodLowerEstimate
            periodLowerEstimate = periodUpperEstimate
            periodUpperEstimate = temp
            
            temp2 = ak
            ak = bk
            bk = temp2

        fractionalErrorOld = fractionalError
        fractionalError = 1.0 * np.abs(period-periodUpperEstimate)/period
        fractionalErrorChange = np.abs(fractionalError - fractionalErrorOld)
        print ak
        print bk 
        print periodLowerEstimate
        print periodUpperEstimate
            
    if fractionalError <= 0.001:
        solution['initialperiod'] = bk
        finalProcessString, finalOutf = buildCommandString(array, solution['initialperiod'][0], 'final')
        
        if finalOutf == -1:
            solution['initialperiod'] = -1
            solution['OBLI'] = -1
            return solution
        
        datafinal = readEvolution(finalOutf)
        maskindicesfinal = np.isclose(datafinal['t'], age)
        solution['OBLI'] = np.mean(datafinal['convincl'][maskindicesfinal])
    else:
        solution['initialperiod'] = -1
        solution['OBLI'] = -1

    return solution


testplanet = np.array([((0.479329, 1.165, 3.405909, 4, 0.0, 8, 5, 1.4, 0.155, 12.0, 2.45, 2.5))], \
       dtype=np.dtype({'names': ('PMASS', 'STARMASS', 'ORBPERIOD', 'AGE', 'OBLI', 'lgQ', 'lgQinr', 'p-disk', 'K', 'tcoup', 'wsat', 'tdisk'), \
        'formats':('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8' )}))
solution = solve(testplanet)
print solution
