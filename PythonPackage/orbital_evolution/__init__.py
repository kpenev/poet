"""The module for calculating tidal orbital evolution."""

from glob import glob
import os.path

from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from orbital_evolution.initial_condition_solver import InitialConditionSolver
from orbital_evolution.transformations import phase_lag, lgQ
from orbital_evolution.binary import Binary
from orbital_evolution.learning_ic_solver import LearningICSolver

__all__ = [os.path.basename(fname).rstrip('.py')
           for fname in glob(os.path.join(__path__[0], '*.py'))]
