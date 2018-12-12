LIBRARIES = $(shell (find ../../build/libs/*/shared/debug -maxdepth 1 -type f -name "*lib*"))
LIB_PATHS = $(dir $(LIBRARIES))
LIB_NAMES = $(patsubst lib%,%, $(basename $(notdir $(LIBRARIES))))

EXECUTABLES = $(basename $(shell (ls *.cpp)))

CPP_FLAGS = -std=c++11 -g -Wall -O0 \
			-DTOOLCHAIN_GCC \
			-I../third_party_libs/Eigen/ \
			-I../third_party_libs/alglib/ \
			-I../third_party_libs/dirent_hacked/ \
			$(addprefix -L, $(LIB_PATHS)) \
			../Core/*.cpp \
			../third_party_libs/alglib/alglib/src/*.cpp
LINK_FLAGS = $(addprefix -l, $(LIB_NAMES)) \
			 -lgsl \
			 -lgslcblas \
			 -lboost_serialization \
			 -lpthread \
			 -lm \
			 -m64

all: $(EXECUTABLES)

%: %.cpp
	g++  $(CPP_FLAGS) -o $@ $@.cpp $(LINK_FLAGS)

clean:
	rm $(EXECUTABLES)
