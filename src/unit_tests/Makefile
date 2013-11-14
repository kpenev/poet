CPP=g++
ALGLIB_OBJ=../alglib/src/interpolation.o ../alglib/src/ap.o \
		   ../alglib/src/alglibinternal.o ../alglib/src/optimization.o \
		   ../alglib/src/linalg.o ../alglib/src/integration.o \
		   ../alglib/src/alglibmisc.o ../alglib/src/solvers.o \
		   ../alglib/src/specialfunctions.o
COMBINED_CODE = ../Functions.o ../OrbitSolver.o ../Planet.o ../Star.o ../StellarEvolution.o ../YRECIO.o ../StellarSystem.o ../Common.o ../StoppingConditions.o Common.o
TEST_LIB = -lcpptest -lgsl -lblas -lboost_serialization-mt
ALL_EXEC=test_InterpolatingFunctionALGLIB test_EvolvingStellarQuantity test_Star test_Planet test_StellarEvolution test_StellarSystem
CPPFLAGS = -g -Wall -I../ -I/opt/local/include -L/opt/local/lib -D STANDALONE 

default: all
Common.o: Common.h ../*.h

%: %.cpp %.h Common.o ../*.cpp ../*.h ../*.o ../alglib/src/*.cpp ../alglib/src/*.h
#	make -C ../
	$(CPP) $(CPPFLAGS) $< -c -o $@.o
	$(CPP) $(CPPFLAGS) $(COMBINED_CODE) $(ALGLIB_OBJ) $@.o -o $@ $(TEST_LIB)

all: $(ALL_EXEC)

clean: 
	rm -f $(ALL_EXEC) *.o 
