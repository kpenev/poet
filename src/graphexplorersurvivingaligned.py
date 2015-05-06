import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas
import scipy
import fnmatch

import os
import sys


#need to cluster planets
logfile = open('./explorationpart/SurvivingAlignedPlanets/Survivingalignedgraphlogs.txt', 'w', 0)

incl = np.linspace(0, np.pi * 0.98, 16)
logQinr = np.array([3, 5])
logQ = np.array([5, 7, 9])
mstar = np.array([1, 0.5])
mplanet = np.array([1])
orbperiods = np.array([3, 4, 5, 6])

fastrot = np.array([(1.4, 0.155, 12, 2.45, 2.5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))
mediumrot = np.array([(7, 0.17, 28, 2.45, 5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))
slowrot = np.array([(10, 0.17, 30, 2.45, 5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))

rot = [fastrot, mediumrot, slowrot]
rottxt = ['fast', 'medium', 'slow']


planetlist = []
ageslist = np.array([0.01*1.5**i for i in range(18)])
ageslist = np.around(ageslist, 5)
allEvolutions = os.listdir('./explorationpart/')


def readEvolution(location):
    location = './explorationpart/' + location
    #output columns will be t, period, a, convincl, radincl, Lconv, Lrad, Iconv, Irad, I, mode
    data = np.genfromtxt(location, skiprows=1, usecols=(0, 1, 2, 3, 4, 10), dtype={'names': ('t', 'ORBPERIOD', 'a', 'convincl', 'radincl', 'mode'), 'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'S15')})
    return data

#rotation should be fast.evo, medium.evol or slow.evol
def buildMatchingString(M, m, period, Qin, Qout, Rotation):
    string = 'M_' + str(M) + '*' + 'm_' + str(m) + '*' + 'period_' + str(period) + '*' + 'lgQin_' + str(Qin) + '*' + \
             'lgQout_' + str(Qout) + '*' + 'Obli_*' + 'rot_' + Rotation + '.evol'
    return string

def buildGraphString(M, m, period, Qin, Qout, Rotation):
    string = 'SurvivingAlignedM_' + str(M) + 'm_' + str(m) + 'period_' + str(period)  + 'lgQin_' + str(Qin)  + \
             'lgQout_' + str(Qout) + 'rot_' + Rotation + '.png'
    return string


#setting all evolutions
#evolutions = []
#for f in allEvolutions:
    #if fnmatch.fnmatch(f, buildMatchingString(1.0, 1, 3, 3, 5, 'fast')):
        #evolution = readEvolution(f)
        #evolution['t'] = np.around(evolution['t'], 5)
        #evolutions.append(evolution)

#print evolutions

#print buildMatchingString(1.0, 1, 3, 3, 5, 'fast')
#evolutions = evolutions[0]
#print evolutions['t']

#for evol in evolutions:
    #for i in range(ageslist.shape[0]):
        #obliquityAtAge = np.mean(evolution['convincl'][ageMask])
        #print obliqui`tyAtAge

def buildFractionGraph(M, m, period, Qin, Qout, Rotation):
    #find evolutions
    Evolutions = []
    for f in allEvolutions:
        if fnmatch.fnmatch(f, buildMatchingString(M, m, period, Qin, Qout, Rotation)):
            evolution = readEvolution(f)
            evolution['t'] = np.around(evolution['t'], 5)
            Evolutions.append(evolution)


    fractionAtAges = np.zeros(len(ageslist))

    for i in range(len(ageslist)):
        deadPlanets = 0
        alignedCount = 0
        for evol in Evolutions:  
            ageMask = np.isclose(evol['t'], ageslist[i])
            total = len(Evolutions)
            obliquityAtAge = -1
            #find out which ones are NA
            whereNA = np.isnan(evol['a'][ageMask])
            if np.all(whereNA):
                obliquityAtAge = np.nan
            else: 
                obliquityAtAge = np.mean(evol['convincl'][ageMask])

            if np.isnan(obliquityAtAge):
                deadPlanets = deadPlanets + 1

            elif obliquityAtAge <= np.radians(10):
                alignedCount = alignedCount + 1
            
            
        if deadPlanets == len(Evolutions):
                fractionAtAges[i] = 0
        else:
            fractionAtAges[i] = 1.0 * alignedCount/(len(Evolutions) - deadPlanets)
    
    fig = plt.figure()
    plt.title(r'Fraction of Surviving Aligned Planets')
    data = pandas.DataFrame({'t': ageslist, 'fraction': fractionAtAges})
    axs = sns.regplot(data['t'], data['fraction'])
    axs.set_xlabel('Ages (Gyr)')
    axs.set_ylabel('Fraction')
    axs.set_xscale('log')
    axs.set_xlim([0, 10])
    axs.set_ylim([0, 1.0])
    sns.despine()
    plt.savefig('./explorationpart/SurvivingAlignedPlanets/' + buildGraphString(M, m, period, Qin, Qout, Rotation))
    plt.close(fig)

for M in mstar:
    logfile.write('M: {0}\n'.format(M))
    for m in mplanet:
        logfile.write('m: {0}\n'.format(m))
        for period in orbperiods:
            logfile.write('Period: {0}\n'.format(period))
            for Qin in logQinr:
                logfile.write('Qin: {0}\n'.format(Qin))
                for Qout in logQ:
                    logfile.write('Qout: {0}\n'.format(Qout))
                    for Rotation in rottxt:
                        logfile.write('Rotation: {0}\n'.format(Rotation))
                        buildFractionGraph(M, m, period, Qin, Qout, Rotation)
