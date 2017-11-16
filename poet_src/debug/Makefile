LIBRARIES = $(shell (find ../build/libs/*/shared/debug -maxdepth 1 -type f -name "*lib*"))
LIB_PATHS = $(dir $(LIBRARIES))
LIB_NAMES = $(patsubst lib%,%, $(basename $(notdir $(LIBRARIES))))

EXECUTABLES = $(basename $(shell (ls *.cpp)))

CPP_FLAGS = -std=c++11 -g -Wall -O0 \
			-I../third_party_libs/Eigen/ \
			-I../third_party_libs/alglib/ \
			$(addprefix -L, $(LIB_PATHS)) \
			$(addprefix -l, $(LIB_NAMES))

all: $(EXECUTABLES)

%: %.cpp
	g++  $(CPP_FLAGS) -o $@ $@.cpp

clean:
	rm $(EXECUTABLES)
