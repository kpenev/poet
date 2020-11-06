import sqlite3
from sqlalchemy import create_engine

import cheb_expand_data_model as model
import calculateCoeffs as pms


def main(brandNew = 0):
	
	db_engine = create_engine('sqlite:///cheb_expansions.db',echo=True)
	
	#if brandNew:
		# Assume we're creating the table
	
	for table in model.DataModelBase.metadata.sorted_tables:
		if not db_engine.has_table(table.name):
			table.create(db_engine)
	
	#conn = sqlite3.connect('cheb_expansions.db')
	#c = conn.cursor()

	pms.main()

	#bob = model.PmsVAcc(m=0,s=1,accur=7)
	
	#conn.close()


if __name__ == '__main__':
    main()