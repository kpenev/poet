from distutils.core import setup, Extension
from subprocess import Popen, PIPE
from os.path import join as join_paths
from glob import glob
import sys

gsl_config=Popen(['gsl-config', '--prefix'], stdout=PIPE)
gsl_prefix=gsl_config.communicate()[0].strip()

orbital_evolution_code=['Common.cpp', 'Functions.cpp', 'OrbitSolver.cpp',
                        'Planet.cpp', 'Star.cpp', 'StellarEvolution.cpp',
                        'StellarSystem.cpp', 'YRECIO.cpp']
alglib_code=['alglib/src/interpolation.cpp', 'alglib/src/ap.cpp',
             'alglib/src/alglibinternal.cpp', 'alglib/src/optimization.cpp',
             'alglib/src/linalg.cpp', 'alglib/src/integration.cpp',
             'alglib/src/alglibmisc.cpp', 'alglib/src/solvers.cpp',
             'alglib/src/specialfunctions.cpp']

YREC_tracks=glob(join_paths('YREC', '*.track'))

setup(name='OrbitSolver',
      version='1.0',
      author='Kaloyan Penev',
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
