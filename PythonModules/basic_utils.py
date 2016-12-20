from sqlalchemy.orm import sessionmaker
from contextlib import contextmanager
from tempfile import mkdtemp
from shutil import rmtree

Session = sessionmaker()

@contextmanager
def db_session_scope() :
    """Provide a transactional scope around a series of operations."""
    session = Session()
    try :
        yield session
        session.commit()
    except :
        session.rollback()
        raise
    finally :
        session.close()


@contextmanager
def tempdir_scope() :
    """Create a temporary directory and destroy when no longer needed."""

    dirname = mkdtemp()
    try :
        yield dirname
    finally :
        rmtree(dirname)

class Structure :
    """An empty class used only to hold user defined attributes."""

    def __init__(self, **initial_attributes) :
        """Create a class with (optionally) initial attributes."""

        for attribute_name, attribute_value in initial_attributes.items() :
            setattr(self, attribute_name, attribute_value)

    def format(self, prefix='') :
        """Generate a tree-like representation of self."""

        result=''
        for attr_name in dir(self) :
            if attr_name[0]!='_' :
                attribute=getattr(self, attr_name)
                if isinstance(attribute, Structure) :
                    result+=(prefix
                             +
                             '|-'
                             +
                             attr_name
                             +
                             '\n'
                             +
                             attribute.format(prefix + '| '))
                else : result+=(prefix
                                +
                                '|-'
                                +
                                attr_name
                                +
                                ': '
                                +
                                str(attribute)
                                +
                                '\n')
        return result
