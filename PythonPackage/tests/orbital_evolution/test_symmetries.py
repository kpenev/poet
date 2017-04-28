#!/usr/bin/python3 -u

"""Test some evolution symmetries (i.e. exchangi primary<->secondary) etc."""

import matplotlib
matplotlib.use('TkAgg')

import os.path
import sys
sys.path.insert(0, os.path.abspath('../../'))

from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.evolve_interface import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
import numpy
from astropy import units, constants

if __name__ == '__main__' :
    pass
