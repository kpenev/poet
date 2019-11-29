************
Installation
************


Only Linux and MacOS is currently supported. Sorry Windows has no clean way
getting the GNU science library.

If you run into trouble with the installation steps below, please contact
`Kaloyan Penev <mailto:kaloyan.penev@utdallas.edu>`_.

Installation Steps:
===================

  1. Install  `JAVA <https://www.java.com/en/download/>`_
  
  2. Install required libraries:
  
     2.1 gls

     2.2 glscblas

     2.3 boost_serialization

     2.4 pthread

  3. Run a terminal and ``cd`` to a place where you would like to download POET

  3. Clone the poet repository: ``git clone git@github.com:kpenev/poet.git``

  4. Change to the root directory: ``cd poet``

  5. Run ``./gradlew build``. This should take a few minutes and end with a
     message that the build was successful. Something along the lines of::

         (base) domat:poet kpenev$ ./gradlew build
         Starting a Gradle Daemon (subsequent builds will be faster)

         BUILD SUCCESSFUL in 4m 58s
         72 actionable tasks: 52 executed, 20 up-to-date

  6. The step above will create 4 shared libraries, with paths relative to poet:
     
     6.1 ``build/libs/evolve/shared/release/libevolve.dylib`` (MacOS) or
     ``build/libs/evolve/shared/release/libevolve.so`` (Linux)

     6.2 ``build/libs/planet/shared/release/libplanet.dylib`` (MacOS) or
     ``build/libs/planet/shared/release/libplanet.so`` (Linux)

     6.3 ``build/libs/star/shared/release/libstar.dylib`` (MacOS) or
     ``build/libs/star/shared/release/libstar.so`` (Linux)

     6.4
     ``build/libs/stellarEvolution/shared/release/libstellarEvolution.dylib``
     (MacOS) or
     ``build/libs/stellarEvolution/shared/release/libstellarEvolution.so``
     (Linux)

     You need to copy (or symlink) these to locations where your operating
     system can find them when requested. If you have sudo privileges, the
     standard is to copy them to ``/usr/local/lib/``. Other options are
     operating system specific.

  7. Copy (or symlink) the sub-directory ``PythonPackage`` to a location
     where the python import mechanism can find it. The examples in this
     documentation will assume that the top level of the package is called
     ``poet``. For example if using python3.7 under anaconda::

         cp -r PythonPackage ~/anaconda3/lib/python3.7/site-packages/poet

     Alternatively, you can copy the python package wherever you want and add
     the location to your PYTHONPATH environment variable. For example::

         mkdir -p ~/my_python_packages
         cp -r PythonPackage  ~/my_python_packages/poet
         export PYTHONPATH=$PYTHONPATH:~/my_python_packages

     You may wish to add that last line to your ``~/.bashrc`` file under Linux
     or ``~/.bash_profile`` file under MacOS, otherwise you will have to execute
     that in each shell where you plan to run POET.

  8. Create the defalut stellar evolution interpolator::

         cd scripts
         ./add_default_interpolator

     Caution: this will take many hours, so please plan accordingly. By default
     this script will use 4 parallel processes, if you have more CPUs or if you
     want to leave some unused so you can do other things on your computer while
     his is going on, you can use the --num-threads argument.
