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

def main(em = 0, es = 1):
	
	# Check if the appropriate file already exists; if not, create it
	for table in model.DataModelBase.metadata.sorted_tables:
		if not db_engine.has_table(table.name):
			table.create(db_engine)

	# Calculate some expansions for a specified coefficient
	chebCoeffs, Accurs = pms.main(m=em,s=es)
	
	# Prep our lists of values to be put in the database
	PmsCoeffs = []
	Expansions = []
	
	# Prepare to figure out where we're starting with ids
	biggest = 0

	# Open a db session to gain access to relevant i/o utilities
	with db_session_scope() as db_session:
		
		# Find the biggest id so far
		for bob in db_session.query(model.MAndSToAccuracy):
			if bob.id > biggest:
				biggest = bob.id
		
		# Loop through the data we generated, prep it for the database
		for i in range(len(Accurs)):
			PmsCoeffs.append(model.MAndSToAccuracy(id=biggest + 1 + i,m=em,s=es,
						accuracy = Accurs[i],timestamp = datetime.now(timezone.utc)))
			for y in range(len(chebCoeffs[i])):
				Expansions.append(model.ChebExpansionCoeffs(p_id = biggest + 1 + i,place_in_expansion=y,coefficient_value = chebCoeffs[i][y],
				timestamp = datetime.now(timezone.utc)))
		
		# Add everything to the database
		db_session.add_all(PmsCoeffs)
		db_session.add_all(Expansions)


if __name__ == '__main__':
    main()