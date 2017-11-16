#!/usr/bin/python3 -u

import sys
sys.path.append('../PythonPackage')
from stellar_evolution import Session
from stellar_evolution.basic_utils import db_session_scope
from stellar_evolution.manager_data_model import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

def define_varchange_dependent_variables(db_session) :
    """
    Define the dependent variables for stellar evolution variable change.

    Args:
        - db_session:
            The currently active database session.

    Returns: None
    """

    db_variables = [
        VarchangeDependentVariable(id = id, name = name)
        for id, name in enumerate(['Teff', 'logg', 'L', 'rho'])
    ]
    db_session.add_all(db_variables)

if __name__ == '__main__' :
    db_engine = create_engine('sqlite:///test.sqlite',
                              echo = True)
    Session.configure(bind = db_engine)

    with db_session_scope() as db_session :
        for table in DataModelBase.metadata.sorted_tables :
            if db_engine.has_table(table.name) : continue
            if table.name.startswith('varchange_') :
                table.create(db_engine)
                if table.name == 'varchange_dependent_variables' :
                    define_varchange_dependent_variables(db_session)
