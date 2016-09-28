# Makefile for the POET unit tests.
#
#############################################################################
# Modules from the ALGLIB library

ALGLIB_TARGET   = alglib
ALGLIB          = interpolation ap alglibinternal optimization linalg \
                  integration alglibmisc solvers specialfunctions
ALGLIB_SRCPATH  = ../alglib/src/
ALGLIB_OBJECTS  = $(addprefix $(ALGLIB_SRCPATH), $(addsuffix .o, $(ALGLIB)))

#############################################################################
# Modules from the production code

POET_TARGET     = poet_debug
POET            = EccentricityExpansionCoefficients DissipatingZone \
                  CombinedStoppingCondition SynchronizedCondition \
                  BreakLockCondition SpinOrbitLockInfo OrbitalExpressions \
                  BinarySystem DissipatingBody SecondaryDeathCondition \
                  ZoneOrientation ExponentialDecayDiffRotBody \
                  SaturatingSkumanichWindBody WindSaturationCondition \
                  DiskPlanetSystem MESAIO StellarEvolution Functions Common

POET_OBJECTS    = $(addprefix ../, $(addsuffix .o, $(POET)))

#############################################################################
# Modules common to all unit tests

COMMON          = Common ConstPhaseLagDissipatingZone LaiExpressions \
                  RandomDiskPlanetSystem
COMMON_OBJECTS  = $(addsuffix .o, $(COMMON)) 

#############################################################################
# All available unit tests

UNIT_TESTS = $(basename $(shell ls test_*.cpp))

#############################################################################
# Libraries

BOOST_LIB		= $(shell (echo 'int main() {}' > \
				    	   /tmp/trivial.cpp; g++ /tmp/trivial.cpp \
						   -lboost_serialization -o /tmp/trivial 2>/dev/null) \
					&& \
					echo "-lboost_serialization" \
					|| \
					echo "-lboost_serialization-mt")
GSL_LIB			=$(shell gsl-config --libs 2>/dev/null \
				   || \
				   echo -lgsl -lgslcblas)

#############################################################################
# compiler
COMPILER        = g++
INCLUDE_PATHS   = .. /opt/local/include
INCLUDE_FLAGS   = $(addprefix -I, $(INCLUDE_PATHS))
LINK_PATHS      = /opt/local/lib
LINK_FLAGS      = $(addprefix -L, $(LINK_PATHS))

LIBRARIES		= -lcpptest $(GSL_LIB) $(BOOST_LIB)
CFLAGS_COMMON   = -g -Wall -O0
CFLAGS          = $(CFLAGS_COMMON)

#############################################################################
# Compile production code object files

.PHONY: $(POET_TARGET)

$(POET_TARGET):
	make -C .. debug

#############################################################################
# Compile ALGLIB object files

.PHONY: $(ALGLIB_TARGET)

$(ALGLIB_TARGET):
	@echo "building the ALGLIB library ..."
	@cd $(ALGLIB_SRCPATH) && make COMPILER="$(COMPILER)" CFLAGS="$(CFLAGS)" \
		ALGLIB_MODULES="$(ALGLIB_MODULES)" 2>/dev/null

#############################################################################
# Compile individual tests

test_%: test_%.cpp test_%.h $(COMMON_OBJECTS) $(POET_TARGET) $(ALGLIB_TARGET)
	$(COMPILER) $(CFLAGS) -D STANDALONE $(INCLUDE_FLAGS) $(LINK_FLAGS) -o $@\
                $@.cpp $(COMMON_OBJECTS) $(POET_OBJECTS) $(ALGLIB_OBJECTS) \
				$(LIBRARIES)

#############################################################################
# clean: delete intermediate and target files

clean: 
	@echo "deleting intermediate and target files ..."
	@rm -f $(COMMON_OBJECTS)
	@rm -f $(UNIT_TESTS)
