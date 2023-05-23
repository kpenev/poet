"""Interface to POET's stellar evolution interpolation."""

from glob import glob
import os.path

from sqlalchemy.orm import sessionmaker
from stellar_evolution.MIST_realtime.build_model import \
    MISTTrackMaker,\
    TemporaryMESAWorkDirectory

#This is actually a class not a variable
#pylint: disable=invalid-name
Session = sessionmaker()
#pylint: enable=invalid-name

__all__ = [os.path.basename(fname).rstrip('.py')
           for fname in glob(os.path.join(__path__[0], '*.py'))]
