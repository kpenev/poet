# Makefile for POET, the other executables and the documentation.
#
#############################################################################
# install dir

BINDIR = ../../bin/

#############################################################################
# sources

# modules
#
COMMON_MODULES  = Common Functions OrbitSolver Planet Star StellarEvolution \
                  StellarSystem YRECIO

LIBRARIES       = -lgsl -lgslcblas -lboost_serialization-mt -largtable2

#POET

POET            = poet

POET_MODULES    = $(POET) $(COMMON_MODULES)
POET_OBJECTS    = $(addsuffix .o, $(POET_MODULES))
POET_DEPS       = $(addsuffix .d, $(POET_MODULES))

#SimOne
SIMONE          = SimOne

SIMONE_MODULES  = $(SIMONE) $(COMMON_MODULES)
SIMONE_OBJECTS  = $(addsuffix .o, $(POET_MODULES))
SIMONE_DEPS     = $(addsuffix .d, $(POET_MODULES))

#list of executables
EXECUTABLES     = $(POET) $(SIMONE)
TARGETS         = $(POET) $(SIMONE)

###############################################################################
# alglib:

ALGLIB          = alglib
ALGLIBSRCPATH   = alglib/src/
ALGLIBMODULES   = alglibinternal alglibmisc ap dataanalysis diffequations \
                  fasttransforms integration interpolation linalg \
                  optimization solvers specialfunctions statistics
ALGLIBOBJS      = $(addsuffix .o, $(ALGLIBMODULES))
ALGLIBOBJECTS   = $(addprefix $(ALGLIBSRCPATH), $(ALGLIBOBJS))

###############################################################################
# Doxygen documentation:

DOXY_CONF        = ../documentation/DoxygenConfig
DOXY_DEST        = ../documentation/doxygen
DOXY_HTML        = $(addsuffix /html, $(DOXY_DEST))
DOXY_WWW         = kpenev@huffy.astro.princeton.edu:~/WWW/public/tidal_orbital_evolution/

#############################################################################
# compiler

COMPILER        = g++

INCLUDES        = 
CFLAGS_COMMON   = -W -Wall $(INCLUDES) -D'GLOBAL_SVN="$(shell svnversion -n)"'

CFLAGS_DEBUG   = -ggdb3 -DDEBUG $(CFLAGS_COMMON)
CFLAGS_RELEASE = -O3 $(CFLAGS_COMMON)
CFLAGS_PROFILE = -g -pg $(CFLAGS_RELEASE)

LFLAGS_COMMON   = $(LIBRARIES)

LFLAGS_RELEASE  = $(CFLAGS_RELEASE) $(LFLAGS_COMMON)
LFLAGS_DEBUG    = $(CFLAGS_DEBUG) $(LFLAGS_COMMON)
LFLAGS_PROFILE  = $(CFLAGS_PROFILE) $(LFLAGS_COMMON)

#############################################################################
# default target: build everything

.PHONY: all

all: override CFLAGS = $(CFLAGS_RELEASE)
all: override LFLAGS = $(LFLAGS_RELEASE)

all: $(EXECUTABLES) ctags

.PHONY: debug

debug: override CFLAGS = $(CFLAGS_DEBUG)
debug: override LFLAGS = $(LFLAGS_DEBUG)

debug: $(EXECUTABLES) ctags

.PHONY: profile

profile: override CFLAGS = $(CFLAGS_PROFILE)
profile: override LFLAGS = $(LFLAGS_PROFILE)

profile: $(EXECUTABLES) ctags

#############################################################################

.PHONY: $(POET)

$(POET): override CFLAGS = $(CFLAGS_RELEASE)
$(POET): override LFLAGS = $(LFLAGS_RELEASE)

$(POET): $(POET_OBJECTS) $(ALGLIB)
	@$(COMPILER) -o $@ $(POET_OBJECTS) $(ALGLIBOBJECTS) $(LFLAGS) \
		|| [ $$? -eq 0 ]
	@echo "$(POET) compiled successfully!"

#############################################################################

.PHONY: $(SIMONE)

$(SIMONE): override CFLAGS = $(CFLAGS_RELEASE)
$(SIMONE): override LFLAGS = $(LFLAGS_RELEASE)

$(SIMONE): $(SIMONE_OBJECTS) $(ALGLIB)
	@$(COMPILER) -o $@ $(SIMONE_OBJECTS) $(ALGLIBOBJECTS) $(LFLAGS) \
		|| [ $$? -eq 0 ]
	@echo "$(SIMONE) compiled successfully!"

#############################################################################
# ALGLIB library

.PHONY: $(ALGLIB)

$(ALGLIB):
	@echo "building the ALGLIB library ..."
	@cd $(ALGLIBSRCPATH) && make

#############################################################################
# clean: delete intermediate and target files

.PHONY: clean

clean:
	@echo "deleting intermediate and target files ..."
	@rm -f $(POET_OBJECTS) $(POET_DEPS)
	@rm -f $(SIMONE_OBJECTS) $(SIMONE_DEPS)
	@rm -f $(EXECUTABLES)

#############################################################################
# del: delete everything

.PHONY: del

del:
	@echo "cleaning up ..."
	@make clean
	@cd $(ALGLIBSRCPATH) && make clean

#############################################################################
# build ctags index

.PHONY: ctags

ctags:
	@echo updating ctags index ...
	@ctags -R --c++-types=+fx --extra=+q --excmd=pattern --exclude=Makefile

#############################################################################
# build source

.cpp.o:
	$(COMPILER) $(CFLAGS) -c -MMD -MP $<
	@sed -i -e '1s,\($*\)\.o[ :]*,\1.o $*.d: ,' $*.d

-include $(POET_DEPS)
-include $(SIMONE_DEPS)

#############################################################################
# Prepare and install documentation to Penev's homepage

.PHONY: doc

doc:
	@doxygen $(DOXY_CONF) > /dev/null \
		|| [ $$? -eq 0 ]
	@echo "Documentation generated successfully!"
	@chmod -R a+r $(DOXY_DEST)
	@chmod a+x $(addsuffix /search, $(DOXY_HTML))
	@rsync -azr $(addsuffix /*, $(DOXY_HTML)) $(DOXY_WWW) \
		|| [ $$? -eq 0 ]
	@echo "Documentation uploaded to homepage successfully!"

#############################################################################
# help screen

.PHONY: help
help:
	@echo "Makefile for the P(lanetary) O(rbital) E(volution) due to"\
          "T(ides) package"
	@echo "External libraries needed: gsl, lapack, blas, argtable2."
	@echo ""
	@echo "Makefile targets :"
	@echo "all              : build all executables (default)"
	@echo "debug            : build all executables in debug mode"
	@echo "profile          : build all executables to be used with gprof"
	@echo "install          : TODO build all and copy executables to bin"
	@echo "poet             : build the $(POET) executable"
	@echo "SimOne           : build the $(SIMONE) executable"
	@echo "clean            : delete intermediate files and executables"
	@echo "del              : same as clean, plus clean extern libs"
	@echo "ctags            : rebuild ctags index file"
	@echo ""
