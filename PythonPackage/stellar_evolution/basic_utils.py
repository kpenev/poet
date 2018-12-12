"""Define DB sesssion and temp directory scopes."""

from contextlib import contextmanager
from tempfile import mkdtemp
from shutil import rmtree
from . import Session

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


@contextmanager
def tempdir_scope():
    """Create a temporary directory and destroy when no longer needed."""

    dirname = mkdtemp()
    try:
        yield dirname
    finally:
        rmtree(dirname)
