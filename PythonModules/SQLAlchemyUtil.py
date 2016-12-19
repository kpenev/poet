from sqlalchemy.orm import sessionmaker
from contextlib import contextmanager

Session = sessionmaker()

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


