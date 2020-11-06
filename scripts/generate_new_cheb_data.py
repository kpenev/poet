import sqlite3
from contextlib import contextmanager
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
#from sqlalchemy.dialects.sqlite import DATETIME
from datetime import datetime, timezone
import numpy as np

import cheb_expand_data_model as model
import calculateCoeffs as pms

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

def main(brandNew = 0):

	#if brandNew:
		# Assume we're creating the table
	
	for table in model.DataModelBase.metadata.sorted_tables:
		if not db_engine.has_table(table.name):
			table.create(db_engine)
	
	#conn = sqlite3.connect('cheb_expansions.db')
	#c = conn.cursor()

	chebCoeffs, Accurs = pms.main()
	
	#print(chebCoeffs)
	#print('\n\n\n hi \n\n\n')
	#print(Accurs)
	
	#Expansions = [
	#		model.ChebCo(chebOrder = y,
	#				value = i[y])
	#		for y in range(len(i))
	#	for i in chebCoeffs
	#]
	Expansions = []
	for i in range(len(chebCoeffs)):
		for y in range(len(chebCoeffs[i])):
			Expansions.append(model.ChebCo(chebOrder=y,value = chebCoeffs[i][y],
			timestamp = datetime.now(timezone.utc)))
	#PmsCoeffs = [
	#	#print(i)
	#	model.PmsVAcc(m=0,s=1,
	#				accur = Accurs[i.astype(np.int)],
	#				chebCoeffs = Expansions[i.astype(np.int)])
	#	for i in Accurs
	#]
	PmsCoeffs = []
	for i in range(len(Accurs)):
		#print(i)
		PmsCoeffs.append(model.PmsVAcc(m=0,s=1,
					accur = Accurs[i],timestamp = datetime.now(timezone.utc)))#,
					#chebCoeffs = Expansions[i]))

	#bob = model.PmsVAcc(m=0,s=1,accur=7)
	
	#conn.close()
	
	with db_session_scope() as db_session:
		db_session.add_all(PmsCoeffs)
		db_session.add_all(Expansions)


if __name__ == '__main__':
    main()