from distutils.core import setup, Extension
from subprocess import Popen, PIPE
from os.path import join as join_paths
from glob import glob
import sys

gsl_config=Popen(['gsl-config', '--prefix'], stdout=PIPE)
gsl_prefix=gsl_config.communicate()[0].strip()

orbital_evolution_code=['Common.cpp', 'Functions.cpp', 'OrbitSolver.cpp',
                        'Planet.cpp', 'Star.cpp', 'StellarEvolution.cpp',
                        'StellarSystem.cpp', 'YRECIO.cpp',
                        'StoppingConditions.cpp',
                        'ExternalStoppingConditions.cpp',
                        'CustomStellarEvolution.cpp']
alglib_code=['alglib/src/interpolation.cpp', 'alglib/src/ap.cpp',
             'alglib/src/alglibinternal.cpp', 'alglib/src/optimization.cpp',
             'alglib/src/linalg.cpp', 'alglib/src/integration.cpp',
             'alglib/src/alglibmisc.cpp', 'alglib/src/solvers.cpp',
             'alglib/src/specialfunctions.cpp']

YREC_tracks=glob(join_paths('YREC', '*.track'))

setup(name='OrbitSolver',
      version='1.0',
      author='Kaloyan Penev',
      author_email='kpenev@gmail.com',
      description='Calculate the orbital evolution of a planet-star system.',
      long_description=
      """ For a system consisting of a star and a single planet in a circular
      orbit aligned with the stellar rotation, calculate the evolution of the
      semimajor axis and the rotation of the star. If the star has a mass
      below 1.1 solar masses, the rotation of the radiative core and the
      convective envelope are evolved separately. The following physics is
      included:
        * The evolution of the semimajor axis of the orbit due to the tidal
          dissipation in the star.
        * The evolution of the angular momentum of the stellar convective
          envelope by the tidal coupling.
        * The transfer of angular momentum between the stellar convective and
          radiative zones.
        * The effect of the stellar evolution on the tidal dissipation
          efficiency, and stellar core and envelope spins.
        * The loss of stellar convective zone angular momentum to a
          magnetically launched wind. """,
      ext_modules=[Extension('OrbitSolver',
                             ['OrbitSolverModule.cpp']+
                             orbital_evolution_code+
                             alglib_code,
                             include_dirs=['.',
                                           join_paths(gsl_prefix, 'include')],
                             libraries=['gsl', 'gslcblas',
                                        'boost_serialization-mt'],
                             define_macros=[('YREC_TRACK_PATH',
                                             join_paths(sys.exec_prefix,
                                                        'YREC'))]
                            )],
      data_files=[('YREC', YREC_tracks)]
     )
