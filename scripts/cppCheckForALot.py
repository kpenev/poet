import argparse
from contextlib import contextmanager
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
import cheb_expand_data_model as model
import numpy as np
import matplotlib.pyplot as plt
import lib_pms as pms
from multiprocessing import Pool
import random

import sys
#sys.path.append('../PythonPackage/orbital_evolution')
#sys.path.append('..')
sys.path.append('/home/vortebo/ctime/poet')
sys.path.append('/home/vortebo/ctime/poet/PythonPackage')
import PythonPackage.orbital_evolution.evolve_interface as evolve

def make_test_set():
    random.seed(1)
    #element=[0,0,0]
    result=[]
    for i in range(23999):
        result.append([random.choice([-2,0,2]),random.randrange(0,400),random.randrange(0,1000)/1000])
#        result.append(element)

    #print(result)
    return result

def with_poet(em, es, ee):
    coeff = evolve.library.coeff_new(b"pms_db.db",0,False)
    result=evolve.library.coeff_operator(coeff,em,es,ee,0,False)
    evolve.library.coeff_delete(coeff)
    return result

def python_drone(a_val):
    return pms.p_MS(a_val[2],a_val[0],a_val[1],3,0,None,0,2e-9)

def with_python(all_vals):
    procNum=16
    resList = np.array(())
    with Pool(procNum) as p:
        resList=p.map(python_drone,all_vals)
    return resList

def main():
    # Get command line arguments
    #parser = argparse.ArgumentParser(description='Check accuracy of pms from the pms table')
    #parser.add_argument('em',type=int,help='What m (-2,0,2)')
    #args = parser.parse_args()
    # Make sure arguments are in appropriate range
    #if np.abs(args.em) != 2 and args.em != 0:
    #    print("NO your m is WRONG")
    #    args.em=0
    
    try:
        # description
#        random.seed()
        testset = make_test_set()

        poetOut=np.zeros(len(testset))

        for i in range(len(testset)):
             poetOut[i] = with_poet(testset[i][0],testset[i][1],testset[i][2])
             #print("lol")

        pythonOut=with_python(testset)
        #print(poetOut)
        #print(pythonOut)
        difference=poetOut-pythonOut
        print(difference[difference>2e-9])
        #print("lol again")
        # Calculate
#        calculant = pms.getCoefficient(args.em,args.es,midPoints,midPoints.size,3,None,0,2e-9) # x as opposed to midPoints
        
        #accur = np.max( np.abs(np.array(calculant)[midPoints<=cutoff] - interpolant[midPoints<=cutoff]) )
        # Report
 #       print("Resulting accuracy is: ",accur)
        #return accur
  #      compare_pms(x,y)
        #plt.semilogy(x[midPoints<=cutoff],theDiff,'.')
        #plt.savefig("figures.png")
    except:
        print('There has been an... incident.')
        raise
    
if __name__ == '__main__':
    main()
