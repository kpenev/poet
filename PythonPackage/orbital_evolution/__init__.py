"""The module for calculating tidal orbital evolution."""

from glob import glob
import os.path

__all__ = [os.path.basename(fname).rstrip('.py')
           for fname in glob(os.path.join(__path__[0], '*.py'))]
