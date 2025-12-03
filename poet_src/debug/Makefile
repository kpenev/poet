LIBRARIES = $(shell (find ../../build/libs/*/static/release -maxdepth 1 -type f -name "*lib*"))
LIB_PATHS = $(dir $(LIBRARIES))
LIB_NAMES = $(patsubst lib%,%, $(basename $(notdir $(LIBRARIES))))

EXECUTABLES = $(basename $(shell (ls *.cpp)))

CPP_FLAGS = -std=c++11 -g -Wall -O3 \
			-DTOOLCHAIN_GCC \
			-DNDEBUG \
			-I../third_party_libs/Eigen/ \
			-I../third_party_libs/alglib/ \
			-I../third_party_libs/dirent_hacked/ \
                        -I/opt/apps/intel19/python3_9/boost/1.72/include \
                        -I/opt/apps/intel19/gsl/2.7/include \
                        -isystem /opt/intel/compilers_and_libraries_2020.1.217/linux/pstl/include \
                        -isystem /opt/intel/compilers_and_libraries_2020.1.217/linux/daal/include \
                        -isystem /opt/intel/compilers_and_libraries_2020.1.217/linux/tbb/include \
                        -isystem /opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/include \
                        -isystem /opt/intel/compilers_and_libraries_2020.1.217/linux/ipp/include \
                        -isystem /opt/apps/gcc/9.4.0/lib/gcc/x86_64-pc-linux-gnu/9.4.0/../../../../include/c++/9.4.0 \
                        -isystem /opt/apps/gcc/9.4.0/lib/gcc/x86_64-pc-linux-gnu/9.4.0/../../../../include/c++/9.4.0/x86_64-pc-linux-gnu \
                        -isystem /opt/apps/gcc/9.4.0/lib/gcc/x86_64-pc-linux-gnu/9.4.0/../../../../include/c++/9.4.0/backward \
                        -isystem /opt/apps/gcc/9.4.0/lib/gcc/x86_64-pc-linux-gnu/9.4.0/include \
                        -isystem /usr/local/include \
                        -isystem /opt/apps/gcc/9.4.0/include \
                        -isystem /opt/apps/gcc/9.4.0/lib/gcc/x86_64-pc-linux-gnu/9.4.0/include-fixed \
                        -isystem /usr/include \
			$(addprefix -L, $(LIB_PATHS)) \
			../Core/*.cpp \
			../third_party_libs/alglib/alglib/src/*.cpp
LINK_FLAGS = $(addprefix -l, $(LIB_NAMES)) \
                         -L/opt/apps/intel19/python3_9/boost/1.72/lib \
                         -L/opt/apps/intel19/gsl/2.7/lib \
			 -lgsl \
			 -lgslcblas \
			 -lboost_serialization \
			 -lpthread \
			 -lsqlite3 \
			 -lm \
			 -m64

all: $(EXECUTABLES)

%: 
	g++  $(CPP_FLAGS) -o $@ $@.cpp $(LINK_FLAGS)

clean:
	rm $(EXECUTABLES)
