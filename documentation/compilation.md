Compilation {#compilation}
===========
Compiling the (P)lanetary (O)rbital (E)volution due to (T)ides package
-----------------------------------------------------------------

External libraries needed: **gsl**, **gslcblas**, **boost_serialization**, **argtable2**.

If the **gsl-config** executable is in the search path, it is used to figure out the compiler options that provide the search paths to header files and libraries.

We have encountered BOOST serialization libraries named either **boost_serialization** or **boost_serialization-mt** you can change the name by providing **BOOST_LIB='boost_serialization'** on the command line. By default **boost_serialization-mt** is used.

Depending on your set-up, the **BOOST/argtable2** header files and/or libraries may not be in your standard search paths. If you see compile time errors to that effect, please use the **CPATH** and **LIBRARY_PATH** environment variables to tell g++ where to find those.

Makefile targets :
------------------
 * all              : build all executables (default)
 * debug            : build all executables in debug mode
 * profile          : build all executables to be used with gprof
 * install          : TODO build all and copy executables to ../bin/. Use BINDIR='<install path>' on the command line to change where the binaries are installed.
 * release          : creates and uploads a .tgz file from which poet can be installed ot a publicly accessible location.
 * poet             : build the poet executable
 * SimOne           : build the SimOne executable
 * MESAIO           : build the debug_mesaio executable
 * clean            : delete intermediate files and executables
 * del              : same as clean, plus clean extern libs and documentation
 * ctags            : rebuild ctags index file

