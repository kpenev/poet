"""Make sure the POET library is in the python module search path."""

import os.path
import sys

poet_root = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.join(poet_root, 'PythonPackage'))
