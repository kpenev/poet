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

orbital_evolution_code=['../Common.cpp', '../Functions.cpp',
                        '../OrbitSolver.cpp', '../Planet.cpp', '../Star.cpp',
                        '../StellarEvolution.cpp', '../StellarSystem.cpp',
                        '../YRECIO.cpp', '../StoppingConditions.cpp',
                        '../ExternalStoppingConditions.cpp',
                        '../CustomStellarEvolution.cpp']
alglib_code=['../alglib/src/interpolation.cpp', '../alglib/src/ap.cpp',
             '../alglib/src/alglibinternal.cpp',
             '../alglib/src/optimization.cpp', '../alglib/src/linalg.cpp',
             '../alglib/src/integration.cpp', '../alglib/src/alglibmisc.cpp',
             '../alglib/src/solvers.cpp',
             '../alglib/src/specialfunctions.cpp']

YREC_tracks=glob(join_paths('../YREC', '*.track'))

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
                             include_dirs=['..',
                                           join_paths(gsl_prefix, 'include')],
                             libraries=['gsl', 'gslcblas',
                                        'boost_serialization-mt'],
                             define_macros=[('YREC_TRACK_PATH',
                                             join_paths(data_path, 'YREC'))]
                            )],
      data_files=[(data_path, YREC_tracks)]
     )

import poet
chmod(poet.__file__,
      stat.S_IRUSR | stat.S_IXUSR | stat.S_IWUSR | stat.S_IRGRP |
      stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
