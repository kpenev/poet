import numpy as np
import scipy
import os
import sys

logfile = open('log.txt', 'w')



#can't exactly be pi or the simulation cancels out
incl = np.linspace(0, np.pi * 0.98, 10)
logQinr = np.array([4, 5, 6])
logQ = np.array([6, 7, 8])



planetdata = np.loadtxt('initialplanetsearch.txt', skiprows=4, delimiter=',', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12), \
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
    initialentry = './poet' + ' '
    inputinfo = '--Mstar=' + str(array['STARMASS'][0]) + ' ' + '--Mplanet=' + str(array['PMASS'][0]) + ' ' \
            + '--init-orbit-period=' + str(period) + ' ' + '--require-ages=' + str(array['AGE'][0]) + ' '\
            + '--init-inclination=' + str(array['OBLI'][0]) + ' ' + '--lgQ=' + str(array['lgQ'][0]) + ' '\
            + '--lgQinr=' + str(array['lgQinr'][0]) + ' ' + '--p-disk=' + str(array['p-disk'][0]) + ' '\
            + '-K ' + str(array['K'][0]) + ' ' + '--core-env-coupling-timescale=' + str(array['tcoup'][0]) + ' '\
            + '--high-mass-wind-sat-w=' + str(array['wsat'][0]) + ' '+ '--t-disk=' + str(array['tdisk'][0]) + ' ' + '--tmax=' + str(1.01 * array['AGE'][0]) + ' '
    
    outputinfo = '--output-columns=t,Porb,a,convincl,radincl,Lconv,Lrad,Iconv,Irad,I,mode '
    location = 'solver' + solverEnd + '.evol'
    locationinfo = '--output=' + location
    return (initialentry + inputinfo + outputinfo + locationinfo, location)
    


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
    periodUpper = array['ORBPERIOD'][0] * 2
    periodLower = array['ORBPERIOD'][0] * 0.60
    period = array['ORBPERIOD'][0]

    upperProcessString, upperOutf = buildCommandString(array, periodUpper, UPPER)
    lowerProcessString, lowerOutf = buildCommandString(array, periodLower, LOWER)

    
    upperGood = os.system(upperProcessString)
    lowerGood = os.system(lowerProcessString)
    #simulation could not be allocated
    if (upperGood == 6 or lowerGood == 6):
        solution['initialperiod'] = -1
        solution['OBLI'] = -1
        return solution

    #wait for simulations to finish before continuing to work with them
    #exit = [p.wait() for p in process1, process2]

    datatempUpper = readEvolution(upperOutf)
    datatempLower = readEvolution(lowerOutf)

    maskindicesupper = datatempUpper['t'] == age
    maskindiceslower = datatempLower['t'] == age


    #didn't get to age
    #these simulations were never finishing
    if (not np.any(maskindicesupper)):
        while not np.any(maskindicesupper):
            solution['initialperiod'] = -1
            solution['OBLI'] = -1
            return solution
                
    #the estimates are the output f(b) and f(a) respectively, but notice that they are offset from zero
    #the array call is needed as sometimes it returns two
    periodUpperEstimate = np.array(datatempUpper['ORBPERIOD'][maskindicesupper])[0]
    
    if (not np.any(maskindiceslower)):
        periodLowerEstimate = 0
    else:
        periodLowerEstimate = np.array(datatempLower['ORBPERIOD'][maskindiceslower])[0]
    #make sure that solution is inside region
    
    if np.isnan(periodLowerEstimate):
        periodLowerEstimate = 0
    
    if np.isnan(periodUpperEstimate):
        while np.isnan(periodUpperEstimate):
            print 'Loop 5'
            periodUpper = periodUpper * 1.4
            upperProcessString, upperOutf = buildCommandString(array, periodUpper, UPPER)
            os.system(upperProcessString)
            datatempUpper = readEvolution(upperOutf)
            maskindicesupper = datatempUpper['t'] == age
            
            if (not np.any(maskindicesupper)):
                solution['initialperiod'] = -1
                solution['OBLI'] = -1
                return solution
            
            periodUpperEstimate = np.array(datatempUpper['ORBPERIOD'][maskindicesupper])[0]

    while (periodUpperEstimate - period) < 0:
        periodUpper = periodUpper * 1.4
        upperProcessString, upperOutf = buildCommandString(array, periodUpper, UPPER)
        os.system(upperProcessString)
        datatempUpper = readEvolution(upperOutf)
        maskindicesupper = datatempUpper['t'] == age
        #solution wasn't bounded for Hat p 17 b
        #when the step size was increased
        if (not np.any(maskindicesupper)):
            solution['initialperiod'] = -1
            solution['OBLI'] = -1
            return solution
        periodUpperEstimate = np.array(datatempUpper['ORBPERIOD'][maskindicesupper])[0]


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
        process = os.system(sProcessString)
        datatempS = readEvolution(sOutf)
        maskindicesS = datatempS['t'] == age
        if (not np.any(maskindicesS)):
            periodEstimateS = 0
        else:
            periodEstimateS = np.array(datatempS['ORBPERIOD'][maskindicesS])[0]
        
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
            
    if fractionalError <= 0.001:
        solution['initialperiod'] = bk
        finalProcessString, finalOutf = buildCommandString(array, solution['initialperiod'][0], 'final')
        process = os.system(finalProcessString)
        datafinal = readEvolution(finalOutf)
        maskindicesfinal = np.array(datafinal['t'] == age)
        solution['OBLI'] = np.array(datafinal['convincl'][maskindicesfinal])[0]
    else:
        solution['initialperiod'] = -1
        solution['OBLI'] = -1

    return solution


#testplanet = np.array([((0.52, 0.93, 3.36, 11, 0.3490658503988659, 6.0, 4.0, 1.4, 0.155, 12.0, 2.45, 2.5))], \
       #dtype=np.dtype({'names': ('PMASS', 'STARMASS', 'ORBPERIOD', 'AGE', 'OBLI', 'lgQ', 'lgQinr', 'p-disk', 'K', 'tcoup', 'wsat', 'tdisk'), \
        #'formats':('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8' )}))
#solution = solve(testplanet)
#print solution

fastrot = np.array([(1.4, 0.155, 12, 2.45, 2.5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))
mediumrot = np.array([(7, 0.17, 28, 2.45, 5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))
slowrot = np.array([(10, 0.17, 30, 2.45, 5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))

rot = [fastrot, mediumrot, slowrot]

#these resulted in good runs
goodConfigurations = []

goodConfigurationsFile = open('part1goodconfig.txt', 'w')

fastrotationsfile = open('fastrotationsgoodconfig.txt', 'w')
mediumrotationsfile = open('mediumrotationsgoodconfig.txt', 'w')
slowrotationsfile = open('slowrotationsgoodconfig.txt', 'w')

logQinr4file = open('logQinr4GoodConfig.txt', 'w')
logQinr5file = open('logQinr5GoodConfig.txt', 'w')
logQinr6file = open('logQinr6GoodConfig.txt', 'w')

logQ6file = open('logQ6GoodConfig.txt', 'w')
logQ7file = open('logQ7GoodConfig.txt', 'w')
logQ8file = open('logQ8GoodConfig.txt', 'w')

#this array will have how many of the simulations of a certain type for each planet resulted in obliquity
#the structure is Fast, Medium, Slow, logQinr3, logQinr4, logQinr5, 
planetObliquityCounts = []
planetObliquityProportions = []

for planetdataarray in planetdata:
    #age = np.random(planetdataarray['AGE'] - planetdataarray['AGEMINUS'], planetdataarray['AGE'] + planetdataarray['AGEPLUS'])
    age = planetdataarray['AGE']
    obliLower = np.radians(np.abs(planetdataarray['OBLI'] - planetdataarray['OBLIMINUS']))
    obliUpper = np.radians(np.abs(planetdataarray['OBLI'] + planetdataarray['OBLIPLUS']))
    totalObliquities = incl.shape[0]
    
    currentPlanetObliquityCounts = np.array([(planetdataarray['NAMES'], 0, 0, 0, 0, 0, 0, 0, 0, 0)], \
            dtype=np.dtype({'names': ('NAMES', 'Fast', 'Medium', 'Slow', 'logQinr4', 'logQinr5', 'logQinr6', 'logQ6', 'logQ7', 'logQ8'), \
                                             'formats' : ('S10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}))
    currentPlanetObliquityProportions = np.array([(planetdataarray['NAMES'], 0, 0, 0, 0, 0, 0, 0, 0, 0)], \
            dtype=np.dtype({'names': ('NAMES', 'Fast', 'Medium', 'Slow', 'logQinr4', 'logQinr5', 'logQinr6', 'logQ6', 'logQ7', 'logQ8'), \
                                             'formats' : ('S10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}))

    for Qin in logQinr:
        for Qout in logQ:
            for rotparam in range(3):
                for obli in incl:
                    planetConfig = np.array([(planetdataarray['NAMES'],planetdataarray['PMASS'], planetdataarray['STARMASS'], planetdataarray['ORBPERIOD'], age, obli, \
                                            Qout, Qin, rot[rotparam]['p-disk'], rot[rotparam]['K'], rot[rotparam]['tcoup'], rot[rotparam]['wsat'], rot[rotparam]['tdisk'])], \
                                            dtype=np.dtype({'names': ('NAMES', 'PMASS', 'STARMASS', 'ORBPERIOD', 'AGE', 'OBLI', 'lgQ', 'lgQinr', 'p-disk', 'K', 'tcoup', 'wsat', 'tdisk'), \
                                            'formats':('S10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8' )}))
                    print 'Running configuration for: ' + str(planetConfig)
                    logfile.write('Running configuration for: ' + str(planetConfig) + '\n')
                    solution = solve(planetConfig)
                    print str(solution)
                    logfile.write(str(solution) + '\n')
                    solutionObli = solution['OBLI']
                    #the second case is needed because of negatives
                    #this means we found a good solution
                    if (solutionObli >= obliLower and solutionObli <= obliUpper) or (solutionObli >= obliUpper and solutionObli <= obliLower) and solution['initialperiod'] != -1:
                        print 'Solution Found'
                        logfile.write('Solution found\n')
                        goodConfigurations.append(planetConfig)
                        goodConfigurationsFile.write(str(planetConfig) + '\n')

                        if rotparam == 0:
                            print 'Wrote 1'
                            logfile.write('Wrote 1\n')
                            fastrotationsfile.write(str(planetConfig) + '\n')
                            currentPlanetObliquityCounts['Fast'] = currentPlanetObliquityCounts['Fast'] + 1
                        if rotparam == 1:
                            logfile.write('Wrote 2\n')
                            print 'Wrote 2'
                            mediumrotationsfile.write(str(planetConfig) + '\n')
                            currentPlanetObliquityCounts['Medium'] = currentPlanetObliquityCounts['Medium'] + 1
                        if rotparam == 2:
                            logfile.write(str(planetConfig) + '\n')
                            print 'Wrote 3'
                            slowrotationsfile.write(planetConfig)
                            currentPlanetObliquityCounts['Slow'] = currentPlanetObliquityCounts['Slow'] + 1
                        if Qin == 4:
                            logfile.write(str(planetConfig) + '\n')
                            print 'Wrote 4'
                            logQinr4file.write(str(planetConfig) + '\n')
                            currentPlanetObliquityCounts['logQinr4'] = currentPlanetObliquityCounts['logQinr4'] + 1
                        if Qin == 5:
                            logfile.write('Wrote 5\n')
                            print 'Wrote 5'
                            logQinr5file.write(str(planetConfig) + '\n')
                            currentPlanetObliquityCounts['logQinr5'] = currentPlanetObliquityCounts['logQinr5'] + 1
                        if Qin == 6:
                            logfile.write('Wrote 6\n')
                            print 'Wrote 6'
                            logQinr6file.write(str(planetConfig) + '\n')
                            currentPlanetObliquityCounts['logQinr6'] = currentPlanetObliquityCounts['logQinr6'] + 1
                        if Qout == 6:
                            logfile.write('Wrote 7\n')
                            print 'Wrote 7'
                            logQ6file.write(str(planetConfig) + '\n')
                            currentPlanetObliquityCounts['logQ6'] = currentPlanetObliquityCounts['logQ6'] + 1
                        if Qout == 7:
                            logfile.write('Wrote 8\n')
                            print 'Wrote 8'
                            logQ7file.write(str(planetConfig) + '\n')
                            currentPlanetObliquityCounts['logQ7'] = currentPlanetObliquityCounts['logQ7'] + 1
                        if Qout == 8:
                            logfile.write('Wrote 9\n')
                            print 'Wrote 9'
                            logQ8file.write(str(planetConfig) + '\n')
                            currentPlanetObliquityCounts['logQ8'] = currentPlanetObliquityCounts['logQ8'] + 1

    planetObliquityCounts.append(currentPlanetObliquityCounts)

    currentPlanetObliquityProportions['Fast'] = 1.0 * currentPlanetObliquityCounts['Fast'] / (incl.shape[0] * logQinr.shape[0] * logQ.shape[0])
    currentPlanetObliquityProportions['Medium'] = 1.0 * currentPlanetObliquityCounts['Medium'] / (incl.shape[0] * logQinr.shape[0] * logQ.shape[0])
    currentPlanetObliquityProportions['Slow'] = 1.0 * currentPlanetObliquityCounts['Slow'] / (incl.shape[0] * logQinr.shape[0] * logQ.shape[0])
    currentPlanetObliquityProportions['logQinr4'] = 1.0 * currentPlanetObliquityCounts['logQinr4'] / (incl.shape[0] * len(rot) * logQ.shape[0])
    currentPlanetObliquityProportions['logQinr5'] = 1.0 * currentPlanetObliquityCounts['logQinr5'] / (incl.shape[0] * len(rot) * logQ.shape[0])
    currentPlanetObliquityProportions['logQinr6'] = 1.0 * currentPlanetObliquityCounts['logQinr6'] / (incl.shape[0] * len(rot) * logQ.shape[0])
    currentPlanetObliquityProportions['logQ6'] = 1.0 * currentPlanetObliquityCounts['logQ6'] / (incl.shape[0] * len(rot) * logQinr.shape[0])
    currentPlanetObliquityProportions['logQ7'] = 1.0 * currentPlanetObliquityCounts['logQ7'] / (incl.shape[0] * len(rot) * logQinr.shape[0])
    currentPlanetObliquityProportions['logQ8'] = 1.0 * currentPlanetObliquityCounts['logQ8'] / (incl.shape[0] * len(rot) * logQinr.shape[0])

    planetObliquityProportions.append(currentPlanetObliquityProportions)
    
    logfile.write(str(currentPlanetObliquityCounts) + '\n')
    logfile.write(str(currentPlanetObliquityProportions) + '\n')

planetObliquityCounts = np.array(planetObliquityCounts)
planetObliquityProportions = np.array(planetObliquityProportions)

np.savetxt('Part1PlanetObliquityCounts.txt', planetObliquityCounts, header = 'Name Fast Medium Slow logQinr4 logQinr5 logQinr6 logQ6 logQ7 logQ8', fmt='%s')
np.savetxt('Part1PlanetObliquityProportions.txt', planetObliquityProportions, header = 'Name Fast Medium Slow logQinr4 logQinr5 logQinr6 logQ6 logQ7 logQ8', fmt='%s')
np.savetxt('part1goodconfigarray.txt', np.array(goodConfigurations), fmt='%s')
