CPPFLAGS= -g -Wall -L/u/mzzhang/bin/gsl/lib -L/u/mzzhang/bin/boost/lib -L/opt/local/lib
CPP=g++
ALL_OBJ=Functions.o OrbitSolver.o Planet.o Star.o StellarEvolution.o \
	YRECIO.o Common.o StellarSystem.o sampleHist.o \
	sampleRotation.o SimOne.o
ALL_DEP=alglib $(ALL_OBJ)
ALGLIB_OBJ=alglib/src/interpolation.o alglib/src/ap.o \
		   alglib/src/alglibinternal.o alglib/src/optimization.o \
		   alglib/src/linalg.o alglib/src/integration.o \
		   alglib/src/alglibmisc.o alglib/src/solvers.o \
		   alglib/src/specialfunctions.o
ALL_COMPILED=$(ALL_OBJ) $(ALGLIB_OBJ)
LIB=-lgsl -lgslcblas -lboost_serialization-mt

default: all
Functions.o: Functions.h
OrbitSolver.o: OrbitSolver.h
Planet.o: Planet.h AstronomicalConstants.h Functions.h
Star.o: Star.h AstronomicalConstants.h Functions.h StellarZone.h
StellarEvolution.o: StellarEvolution.h Error.h StellarZone.h Functions.h
YRECIO.o: YRECIO.h StellarEvolution.h Error.h
Common.o: Common.h
StellarSystem.o: StellarSystem.h OrbitSolver.cpp
Simulator.o: *.cpp *.h
ConstrainQ.o: *.cpp *.h
sampleHist.o: sampleHist.h
sampleRotation.o: sampleRotation.h
FastRotators.o: *.cpp *.h
SimOne.o: *.cpp *.h

alglib:
	make -C alglib

all: $(ALL_DEP)
	make -C alglib
	$(CPP) $(CPPFLAGS) -o main $(ALL_COMPILED) $(LIB)

doc:
	doxygen documentation/DoxygenConfig
	chmod -R a+r documentation/doxygen
	chmod a+x documentation/doxygen/html/search
	rsync -azr documentation/doxygen/html/* kpenev@huffy.astro.princeton.edu:~/WWW/public/tidal_orbital_evolution/

clean:
	rm -f $(ALL_OBJ) main
	make -C alglib clean
	make -C unit_tests clean
