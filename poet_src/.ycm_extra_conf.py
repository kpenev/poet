"""Configuration so that YCM can function properly."""

import os
import os.path
from glob import glob
from subprocess import Popen, PIPE
import platform

def get_system_paths():
    """Return a list of the directories to search for <...> includes."""

    cpp_compiler = ('g++' if platform.system() == 'Linux' else 'clang')

    DEVNULL = open(os.devnull, 'w')

    compiler_report = Popen(
        [cpp_compiler, '-E', '-v', '-x', 'c++', '-'],
        stdin=DEVNULL,
        stderr=PIPE
    ).communicate()[1].decode().split('\n')

    sys_paths_start = compiler_report.index(
        '#include <...> search starts here:'
    ) + 1

    num_sys_paths = compiler_report[sys_paths_start - 1:].index(
        'End of search list.'
    ) - 1

    DEVNULL.close()

    return [
        path.strip()
        for path in compiler_report[sys_paths_start : sys_paths_start + num_sys_paths]
    ]


#Function name and signature defined by YCM, so not negotiable.
#pylint: disable=invalid-name
#pylint: disable=unused-argument
def Settings(**kwargs):
    """Return the compiler flags YCM should use."""

    flags = ['-x',
             'c++',
             '-Wall',
             '-Wextra',
             '-std=c++11',
             '-fvisibility=hidden',
             '-DTOOLCHAIN_GCC']


    for dirname in get_system_paths():
        flags.extend(['-isystem', dirname])

    for dirname in glob(os.path.join(os.path.dirname(__file__), '*')):
        if os.path.isdir(dirname):
            flags.extend(['-I', os.path.abspath(dirname)])

    for dirname in glob(os.path.join(os.path.dirname(__file__),
                                     'third_party_libs',
                                     '*')):
        if os.path.isdir(dirname):
            flags.extend(['-I', os.path.abspath(dirname)])

    flags.extend(['-L/usr/lib/x86_64-linux-gnu',
                  '-lgsl',
                  '-lgslcblas',
                  '-lm',
                  '-lboost_serialization',
                  '-lpthread',
                  '-m64'])

    return dict(flags=flags)
#pylint: enable=invalid-name
#pylint: enable=unused-argument

if __name__ == '__main__':
    print('\n\t'.join(Settings()['flags']))
