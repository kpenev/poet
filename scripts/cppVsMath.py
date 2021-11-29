import argparse
from contextlib import contextmanager
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
import cheb_expand_data_model as model
import numpy as np
#import matplotlib.pyplot as plt
import lib_pms as pms

import sys
#sys.path.append('../PythonPackage/orbital_evolution')
#sys.path.append('..')
sys.path.append('/home/vortebo/ctime/poet')
sys.path.append('/home/vortebo/ctime/poet/PythonPackage')
import PythonPackage.orbital_evolution.evolve_interface as evolve

eid=0
mine=0
steps=0
cutoff=0

Session = sessionmaker()
db_engine = create_engine('sqlite:///pms_db.db')
Session.configure(bind=db_engine)
@contextmanager
def db_session_scope():
    """Provide a transactional scope around a series of operations."""
    session = Session()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.close()

def load_pms(em, es):
    global eid,mine,steps,cutoff
    coeff = evolve.library.coeff_new(b"pms_db.db",0,False)
    with db_session_scope() as db_session:
        for row in db_session.query(model.Interpolations).filter_by(m=str(em)).filter_by(s=str(es)):
            eid=row.id
            mine=row.min_interp_e
            steps=row.number_of_steps
            cutoff=row.max_checked_e
        
        if 1.0-mine != 0:
            stepsize=(1.0-mine)/(steps-1)
        else:
            raise Exception('Coefficient is always zero, we got this')
        
        x=np.linspace(mine,1,steps)#np.zeros(steps)
        y=np.zeros(steps)
    #okay now like we step through x and use the c++ stuff
    for i in range(x.size):
        y[i]=evolve.library.coeff_operator(coeff,em,es,x[i],0,False)
    
    evolve.library.coeff_delete(coeff)
    
    return x,y

def main():
    global eid,mine,steps,cutoff
    
    # Check if the appropriate file already exists; if not, complain
    for table in model.DataModelBase.metadata.sorted_tables:
        if not db_engine.has_table(table.name):
            print("No table!")
            return 0
    # Get command line arguments
    parser = argparse.ArgumentParser(description='Check accuracy of pms from the pms table')
    parser.add_argument('em',type=int,help='What m (-2,0,2)')
    parser.add_argument('es',type=int,help='What s (0 - 400)')
    args = parser.parse_args()
    # Make sure arguments are in appropriate range
    if np.abs(args.em) != 2 and args.em != 0:
        print("NO your m is WRONG")
        args.em=0
    if args.es < 0 or args.es > 400:
        print("NO your s is WRONG")
        args.es=0
    
    try:
        # Load specified coefficient from database
        x,y = load_pms(args.em,args.es)
        # Find the points we're checking
        newSteps = steps * 2 - 1
        midPoints = np.linspace(mine,1,steps) #[1::2] #newSteps)[1::2]
        # Interpolate
        interpolant = np.interp(midPoints,x,y)
        # Calculate
        calculant = pms.getCoefficient(args.em,args.es,x,x.size,3,None,0,2e-9) # x as opposed to midPoints
        # Compare
        #print(interpolant[7])
        #print(calculant[7])
        #print(interpolant.size)
        print(np.array(calculant).size)
        print(y.size)
        #print( y.tolist())#np.abs(np.array(calculant)[midPoints<=cutoff] - interpolant[midPoints<=cutoff]) )
        accur = np.max( np.abs(np.array(calculant)[midPoints<=cutoff] - y[x<=cutoff]) )
        # Report
        print("Resulting accuracy is: ",accur)
        return accur
    except:
        print('There has been an... incident.')
        raise
    
if __name__ == '__main__':
    main()
