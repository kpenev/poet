CPP=g++
ALGLIB_OBJ=../alglib/src/interpolation.o ../alglib/src/ap.o \
		   ../alglib/src/alglibinternal.o ../alglib/src/optimization.o \
		   ../alglib/src/linalg.o ../alglib/src/integration.o \
		   ../alglib/src/alglibmisc.o ../alglib/src/solvers.o \
		   ../alglib/src/specialfunctions.o
COMBINED_CODE = ../EccentricityExpansionCoefficients.o ../DissipatingZone.o \
				../CombinedStoppingCondition.o ../SynchronizedCondition.o \
				../BreakLockCondition.o ../SpinOrbitLockInfo.o \
				../OrbitalExpressions.o ../BinarySystem.o \
				../DissipatingBody.o ../SecondaryDeathCondition.o \
				../ZoneOrientation.o ../ExponentialDecayDiffRotBody.o \
				../SaturatingSkumanichWindBody.o \
				../WindSaturationCondition.o ../DiskPlanetSystem.o Common.o \
				ConstPhaseLagDissipatingZone.o LaiExpressions.o \
				RandomDiskPlanetSystem.o
TEST_LIB = -lcpptest -lgsl -lblas -lboost_serialization-mt
ALL_EXEC=test_InterpolatingFunctionALGLIB test_EvolvingStellarQuantity \
		 test_Star test_Planet test_StellarEvolution test_StellarSystem \
		 test_EccentricityExpansionCoefficients test_DissipatingZone \
		 test_DissipatingBody test_BinarySystem
CPPFLAGS = -g -Wall -O0 -I../ -I/opt/local/include -L/opt/local/lib -D STANDALONE 

default: all
	
Common.o: Common.h ../*.h

LaiExpressions.o: LaiExpressions.h

ConstPhaseLagDissipatingZone.o: ConstPhaseLagDissipatingZone.h ConstPhaseLagDissipatingZone.cpp Common.h ../*.h

RandomDiskPlanetSystem.o: RandomDiskPlanetSystem.h RandomDiskPlanetSystem.cpp ../*.h

%: %.cpp %.h $(COMBINED_CODE) 
	$(CPP) $(CPPFLAGS) $< -c -o $@.o
	$(CPP) $(CPPFLAGS) $(COMBINED_CODE) $(ALGLIB_OBJ) $@.o -o $@ $(TEST_LIB)

all: $(ALL_EXEC)

clean: 
	rm -f $(ALL_EXEC) *.o 
