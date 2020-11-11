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
	#print(len(Accurs))
	
	biggest = 0
	
	with db_session_scope() as db_session:
		
		# query biggest id so far
		#for id in db_session.query(model.PmsVAcc.id):
		for bob in db_session.query(model.PmsVAcc):
			print(bob.id)
			if bob.id > biggest:
				biggest = bob.id
		
		PmsCoeffs = []
		#for i in range(len(Accurs)):
		#	#print(i)
		#	PmsCoeffs.append(model.PmsVAcc(id=biggest + 1 + i,m=0,s=1,
		#				accur = Accurs[i],timestamp = datetime.now(timezone.utc)))#,
		#				#chebCoeffs = Expansions[i]))
		
		#Expansions = [
		#		model.ChebCo(chebOrder = y,
		#				value = i[y])
		#		for y in range(len(i))
		#	for i in chebCoeffs
		#]
		Expansions = []
		# bob = zip  (pms & chebco)
		# (i,j) in that
		# = pms[i],chebco[i,:]
		#for i,bethune in see above(chebCoeffs):
		#	for index,val in enumerate(chebCoeffs[i]):
		#		Expansions.append(model.ChebCo(p_id = grab the id from the appropriate pmscough,chebOrder=y,value = chebCoeffs[i][y],
		#		timestamp = datetime.now(timezone.utc)))
				
				
		for i in range(len(Accurs)):
			PmsCoeffs.append(model.PmsVAcc(id=biggest + 1 + i,m=0,s=1,
						accur = Accurs[i],timestamp = datetime.now(timezone.utc)))
			for y in range(len(chebCoeffs[i])):
				Expansions.append(model.ChebCo(p_id = biggest + 1 + i,chebOrder=y,value = chebCoeffs[i][y],
				timestamp = datetime.now(timezone.utc)))
		#PmsCoeffs = [
		#	#print(i)
		#	model.PmsVAcc(m=0,s=1,
		#				accur = Accurs[i.astype(np.int)],
		#				chebCoeffs = Expansions[i.astype(np.int)])
		#	for i in Accurs
		#]


		#bob = model.PmsVAcc(m=0,s=1,accur=7)
		
		#conn.close()
		
		
		db_session.add_all(PmsCoeffs)
		db_session.add_all(Expansions)


if __name__ == '__main__':
    main()