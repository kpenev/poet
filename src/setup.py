#!/usr/bin/python

from argparse import ArgumentParser
from distutils.core import setup, Extension
from subprocess import Popen, PIPE
from os.path import join as join_paths
from os import chmod
from glob import glob
import sys
import stat

gsl_config=Popen(['gsl-config', '--prefix'], stdout=PIPE)
gsl_prefix=gsl_config.communicate()[0].strip()

orbital_evolution_code=['src/Common.cpp', 'src/Functions.cpp',
                        'src/OrbitSolver.cpp', 'src/Planet.cpp',
                        'src/Star.cpp',
                        'src/StellarEvolution.cpp', 'src/StellarSystem.cpp',
                        'src/YRECIO.cpp', 'src/StoppingConditions.cpp',
                        'src/ExternalStoppingConditions.cpp',
                        'src/CustomStellarEvolution.cpp']
alglib_code=['src/alglib/src/interpolation.cpp', 'src/alglib/src/ap.cpp',
             'src/alglib/src/alglibinternal.cpp',
             'src/alglib/src/optimization.cpp', 'src/alglib/src/linalg.cpp',
             'src/alglib/src/integration.cpp',
             'src/alglib/src/alglibmisc.cpp',
             'src/alglib/src/solvers.cpp',
             'src/alglib/src/specialfunctions.cpp']

YREC_tracks=glob(join_paths('src/YREC', '*.track'))

data_path="@DATA_PATH@"

setup(name='POET',
      version='0.1.0',
      author='Kaloyan Penev',
      author_email='kpenev@gmail.com',
      description='Calculate the orbital evolution of a planet-star system.',
      long_description=open('README.txt').read(),
      ext_modules=[Extension('poet',
                             ['poetModule.cpp']+
                             orbital_evolution_code+
                             alglib_code,
                             include_dirs=['src',
                                           join_paths(gsl_prefix, 'include')],
                             libraries=['gsl', 'gslcblas',
                                        'boost_serialization-mt'],
                             define_macros=[('YREC_TRACK_PATH',
                                             join_paths(data_path, 'YREC'))],
                             depends=(glob('src/*.h')+
                                      glob('src/alglib/src/*.h'))
                            )],
      data_files=[(join_paths(data_path, 'YREC'), YREC_tracks)]
     )

#import poet
#chmod(poet.__file__,
#      stat.S_IRUSR | stat.S_IXUSR | stat.S_IWUSR | stat.S_IRGRP |
#      stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
