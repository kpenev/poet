import numpy as np
import scipy
import os

def buildCommandString(array, period, solverEnd):
    initialentry = './poet' + ' '
    inputinfo = '--Mstar=' + str(array['STARMASS'][0]) + ' ' + '--Mplanet=' + str(array['PMASS'][0]) + ' ' \
            + '--init-orbit-period=' + str(period) + ' ' + '--require-ages=' + str(array['AGE'][0]) + ' '\
            + '--init-inclination=' + str(array['OBLI'][0]) + ' ' + '--lgQ=' + str(array['lgQ'][0]) + ' '\
            + '--lgQinr=' + str(array['lgQinr'][0]) + ' ' + '--p-disk=' + str(array['p-disk'][0]) + ' '\
            + '-K ' + str(array['K'][0]) + ' ' + '--core-env-coupling-timescale=' + str(array['tcoup'][0]) + ' '\
            + '--high-mass-wind-sat-w=' + str(array['wsat'][0]) + ' '+ '--t-disk=' + str(array['tdisk'][0]) + ' '
    
    outputinfo = '--output-columns=t,Porb,a,convincl,radincl,Lconv,Lrad,Iconv,Irad,I,mode '
    location = 'solver' + solverEnd + '.evol'
    locationinfo = '--output=' + location
    return (initialentry + inputinfo + outputinfo + locationinfo, location)

incl = np.linspace(0, np.pi * 0.98, 16)
logQinr = np.array([5])
logQ = np.array([5, 7, 9])
mstar = np.array([1, 0.5])
mplanet = np.array([1])
orbperiods = np.array([1, 2, 3, 4, 5, 6])

fastrot = np.array([(1.4, 0.155, 12, 2.45, 2.5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))
mediumrot = np.array([(7, 0.17, 28, 2.45, 5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))
slowrot = np.array([(10, 0.17, 30, 2.45, 5)], dtype=np.dtype([('p-disk', 'f8'), ('K', 'f8'), ('tcoup', 'f8'), ('wsat', 'f8'), ('tdisk', 'f8')]))

rot = [fastrot, mediumrot, slowrot]
rottxt = ['fast.evol', 'medium.evol', 'slow.evol']


planetlist = []

for M in mstar:
    for m in mplanet:
        for period in orbperiods:
            for Qin in logQinr:
                for Qout in logQ:
                    for obliquity in incl:
                        for rotationsetting in range(3):
                            outf = './explorationpart/M_' + str(M) + 'm_' + str(m) + 'period_' + str(period) + 'lgQin_' + str(Qin) + \
                                    'lgQout_' + str(Qout) + 'Obli_' + str(obliquity) + 'rot_' + rottxt[rotationsetting]
                            planet = np.array([M, m, period, Qin, Qout, obliquity, rot[rotationsetting]['p-disk'][0], rot[rotationsetting]['K'][0], \
                                    rot[rotationsetting]['tcoup'][0], rot[rotationsetting]['wsat'][0], rot[rotationsetting]['tdisk'][0], outf])
                            planetlist.append(planet)


planetlist = np.array(planetlist)

np.savetxt('./explorationpart/part2exploration.txt', planetlist, header='#M m pform lgQinr lgQ incl pdisk K tcoup wsat tdisk', delimiter=' ', fmt='%s')


outputinfo = '--output-columns=t,Porb,a,convincl,radincl,Lconv,Lrad,Iconv,Irad,I,mode '


beginningofcommand = './poet --input=./explorationpart/part2exploration.txt --input-columns=M,m,pform,lgQinr,lgQ,incl,pdisk,K,tcoup,wsat,tdisk ' 

agesstring = ','.join([str(0.01*1.5**i) for i in range(18)])

os.system(beginningofcommand + '--require-ages=' + agesstring + ' ' + outputinfo)

