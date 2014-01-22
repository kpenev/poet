# Makefile for POET, the other executables and the documentation.
#
#############################################################################
# install dir

BINDIR = ../bin/

#############################################################################
# sources

# modules
#
COMMON_MODULES  = Common Functions OrbitSolver Planet Star StellarEvolution \
                  StellarSystem YRECIO MESAIO StoppingConditions \
				  ExternalStoppingConditions CustomStellarEvolution

BOOST_LIB		= boost_serialization-mt
LIBRARIES       = $(shell gsl-config --libs 2>/dev/null || \
					echo -lgsl -lgslcblas) -l$(BOOST_LIB) -largtable2

#POET

POET            = poet
COMPILE_POET	= compile_poet

POET_MODULES    = $(POET) $(COMMON_MODULES)
POET_OBJECTS    = $(addsuffix .o, $(POET_MODULES))
POET_DEPS       = $(addsuffix .d, $(POET_MODULES)) \
				  $(addsuffix .d-e, $(POET_MODULES))

#SimOne
SIMONE          = SimOne
COMPILE_SIMONE	= compile_SimOne

SIMONE_MODULES  = $(SIMONE) $(COMMON_MODULES)
SIMONE_OBJECTS  = $(addsuffix .o, $(SIMONE_MODULES))
SIMONE_DEPS     = $(addsuffix .d, $(SIMONE_MODULES)) \
				  $(addsuffix .d-e, $(SIMONE_MODULES))

#MESAIO
MESAIO          = debug_mesaio
COMPILE_MESAIO	= compile_MESAIO

MESAIO_MODULES  = $(MESAIO) $(COMMON_MODULES)
MESAIO_OBJECTS  = $(addsuffix .o, $(MESAIO_MODULES))
MESAIO_DEPS     = $(addsuffix .d, $(MESAIO_MODULES)) \
				  $(addsuffix .d-e, $(MESAIO_MODULES))

#list of executables
EXECUTABLES     = $(POET) $(MESAIO) $(SIMONE)
COMPILES        = $(COMPILE_POET) $(COMPILE_MESAIO) $(COMPILE_SIMONE)

###############################################################################
# alglib:

ALGLIB          = alglib
ALGLIB_SRCPATH  = alglib/src/
ALGLIB_MODULES  = alglibinternal alglibmisc ap dataanalysis diffequations \
                  fasttransforms integration interpolation linalg \
                  optimization solvers specialfunctions statistics
ALGLIB_OBJS     = $(addsuffix .o, $(ALGLIB_MODULES))
ALGLIB_OBJECTS  = $(addprefix $(ALGLIB_SRCPATH), $(ALGLIB_OBJS))

###############################################################################
# Doxygen documentation:

DOC_ROOT         = ../documentation
DOXY_CONF        = $(DOC_ROOT)/DoxygenConfig
DOXY_DEST        = $(DOC_ROOT)/doxygen
DOXY_HTML        = $(addsuffix /html, $(DOXY_DEST))
DOXY_WWW         = kpenev@huffy.astro.princeton.edu:~/WWW/public/tidal_orbital_evolution/

#############################################################################
# compiler

COMPILER        = g++

INCLUDES        = 
DATADIR			= .
CFLAGS_COMMON   = -W -Wall $(INCLUDES) \
				  $(shell gsl-config --cflags 2>/dev/null) \
				  -D'GLOBAL_SVN="$(shell svnversion -n)"' \
				  -D'DATADIR="$(DATADIR)"'

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

all: $(COMPILES) ctags

.PHONY: debug

debug: override CFLAGS = $(CFLAGS_DEBUG)
debug: override LFLAGS = $(LFLAGS_DEBUG)

debug: $(COMPILES) ctags

.PHONY: profile

profile: override CFLAGS = $(CFLAGS_PROFILE)
profile: override LFLAGS = $(LFLAGS_PROFILE)

profile: $(COMPILES) ctags

#############################################################################

.PHONY: $(POET)

$(POET): override CFLAGS = $(CFLAGS_RELEASE)
$(POET): override LFLAGS = $(LFLAGS_RELEASE)
$(POET): $(COMPILE_POET)

$(COMPILE_POET): $(POET_OBJECTS) $(ALGLIB)
	@$(COMPILER) -o $(POET) $(POET_OBJECTS) $(ALGLIB_OBJECTS) $(LFLAGS) \
		|| [ $$? -eq 0 ]
	@echo "$(POET) compiled successfully!"

#############################################################################

.PHONY: $(SIMONE)

$(SIMONE): override CFLAGS = $(CFLAGS_RELEASE)
$(SIMONE): override LFLAGS = $(LFLAGS_RELEASE)
$(SIMONE): $(COMPILE_SIMONE)

$(COMPILE_SIMONE): $(SIMONE_OBJECTS) $(ALGLIB)
	@$(COMPILER) -o $(SIMONE) $(SIMONE_OBJECTS) $(ALGLIB_OBJECTS) $(LFLAGS) \
		|| [ $$? -eq 0 ]
	@echo "$(SIMONE) compiled successfully!"

#############################################################################

.PHONY: $(MESAIO)

$(MESAIO): override CFLAGS = $(CFLAGS_RELEASE)
$(MESAIO): override LFLAGS = $(LFLAGS_RELEASE)
$(MESAIO): $(COMPILE_MESAIO)

$(COMPILE_MESAIO): $(MESAIO_OBJECTS) $(ALGLIB)
	@$(COMPILER) -o $(MESAIO) $(MESAIO_OBJECTS) $(ALGLIB_OBJECTS) $(LFLAGS) \
		|| [ $$? -eq 0 ]
	@echo "$(MESAIO) compiled successfully!"

#############################################################################
# ALGLIB library

.PHONY: $(ALGLIB)

$(ALGLIB):
	@echo "building the ALGLIB library ..."
	@cd $(ALGLIB_SRCPATH) && make COMPILER="$(COMPILER)" CFLAGS="$(CFLAGS)" \
		ALGLIB_MODULES="$(ALGLIB_MODULES)" 2>/dev/null

#############################################################################
# deploy executables

.PHONY: install

install: $(EXECUTABLES)
	@echo "installing executables to the bin directory ..."
	@mkdir -p $(BINDIR)
	@mkdir -p $(DATADIR)/YREC
	@sed -e "s%@DATA_PATH@%$(DATADIR)%" setup.py > POETModule/setup.py
	@cp -f $(EXECUTABLES) $(BINDIR)
	@cp YREC/*.track $(DATADIR)/YREC/
	@cp -r serialized_evolution $(DATADIR)
	@cd POETModule && python setup.py install
	cd ..

#############################################################################
# clean: delete intermediate and target files

.PHONY: clean

clean:
	@echo "deleting intermediate and target files ..."
	@rm -f $(POET_OBJECTS) $(POET_DEPS)
	@rm -f $(SIMONE_OBJECTS) $(SIMONE_DEPS)
	@rm -f $(MESAIO_OBJECTS) $(MESAIO_DEPS)
	@rm -f $(EXECUTABLES)
	@cd POETModule && python setup.py clean
	@rm -rf POETModule/build

#############################################################################
# del: delete everything

.PHONY: del

del:
	@echo "cleaning up ..."
	@make clean
	@cd $(ALGLIB_SRCPATH) && make clean ALGLIB_MODULES="$(ALGLIB_MODULES)"
	@rm -rf $(DOXY_HTML)

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
-include $(MESAIO_DEPS)

#############################################################################
# Prepare and install documentation to Penev's homepage

.PHONY: doc

doc: CFLAGS_RELEASE+=-D'COLUMN_NAME_EMPHASIS="__"'

doc: clean $(POET)
	@echo "Compilation {#compilation}" > $(DOC_ROOT)/compilation.md
	@echo "===========" >> $(DOC_ROOT)/compilation.md
	@make help >> $(DOC_ROOT)/compilation.md
	@echo "Usage {#usage}" > $(DOC_ROOT)/usage.md
	@echo "=====" >> $(DOC_ROOT)/usage.md
	@./poet --doxygen-help |sed -f format_doxyhelp.sed >>$(DOC_ROOT)/usage.md
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
	@echo "Compiling the (P)lanetary (O)rbital (E)volution due to (T)ides"\
		  "package"
	@echo "-----------------------------------------------------------------"
	@echo ""
	@echo "External libraries needed: **gsl**, **gslcblas**,"\
		  "**boost_serialization**, **argtable2**."
	@echo ""
	@echo "If the **gsl-config** executable is in the search path, it is"\
		  "used to figure out the compiler options that provide the search"\
		  "paths to header files and libraries."
	@echo ""
	@echo "We have encountered BOOST serialization libraries named either"\
          "**boost_serialization** or **boost_serialization-mt** you can"\
		  "change the name by providing **BOOST_LIB='boost_serialization'**"\
		  "on the command line. By default **boost_serialization-mt** is"\
		  "used."
	@echo ""
	@echo "Depending on your set-up, the **BOOST/argtable2** header files"\
		  "and/or libraries may not be in your standard search paths. If"\
		  "you see compile time errors to that effect, please use the"\
		  "**CPATH** and **LIBRARY_PATH** environment variables to tell g++"\
		  "where to find those."
	@echo ""
	@echo "Makefile targets :"
	@echo "------------------"
	@echo " * all              : build all executables (default)"
	@echo " * debug            : build all executables in debug mode"
	@echo " * profile          : build all executables to be used with gprof"
	@echo " * install          : TODO build all and copy executables to"\
	                             "$(BINDIR). Use BINDIR='<install path>' on"\
							     "the command line to change where the"\
							     "binaries are installed." 
	@echo " * poet             : build the $(POET) executable"
	@echo " * SimOne           : build the $(SIMONE) executable"
	@echo " * MESAIO           : build the $(MESAIO) executable"
	@echo " * clean            : delete intermediate files and executables"
	@echo " * del              : same as clean, plus clean extern libs and"\
							     "documentation"
	@echo " * ctags            : rebuild ctags index file"
	@echo ""
