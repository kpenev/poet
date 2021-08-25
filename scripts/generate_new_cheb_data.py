# Import major libraries
import sqlite3
from contextlib import contextmanager
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from datetime import datetime, timezone

# Import our own python files
import cheb_expand_data_model as model
import calculateCoeffs as pms

# Configure our SQL session
Session = sessionmaker()
db_engine = create_engine('sqlite:///cheb_expansions.db',echo=True)
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

def main(em = 0, es = 1, accGoal = 2e-9):
    
    # Check if the appropriate file already exists; if not, create it
    for table in model.DataModelBase.metadata.sorted_tables:
        if not db_engine.has_table(table.name):
            table.create(db_engine)

    # Calculate some expansions for a specified coefficient
    minE,stepNum,maxE,chebVals = pms.main(m=em,s=es,tolerance=accGoal,kind=3)
    
    # Prep our lists of values to be put in the database
    PmsCoeffs = []
    Expansions = []
    
    # Prepare to figure out where we're starting with ids
    biggest = 0

    # Open a db session to gain access to relevant i/o utilities
    with db_session_scope() as db_session:
        
        # Find the biggest id so far
        for i in db_session.query(model.Interpolations):
            if i.id > biggest:
                biggest = i.id
        
        # Loop through the data we generated, prep it for the database
        PmsCoeffs.append(model.Interpolations(id=biggest + 1,m=em,s=es,min_interp_e=minE,number_of_steps=stepNum,max_checked_e=maxE,
                    interp_accuracy = accGoal,timestamp = datetime.now(timezone.utc)))
        for y in range(len(chebVals)):
            Expansions.append(model.InterpolationData(p_id = biggest + 1,step_number=y,y_value = chebVals[y],
            timestamp = datetime.now(timezone.utc)))
        
        # Add everything to the database
        db_session.add_all(PmsCoeffs)
        db_session.add_all(Expansions)


if __name__ == '__main__':
    main()