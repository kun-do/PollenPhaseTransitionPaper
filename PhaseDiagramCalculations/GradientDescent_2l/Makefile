CXX ?= g++
CXXFLAGS ?= -O3 -std=c++11
EIGEN_DIR ?=
COMMON_DIR = ../Common
SOURCES = main.cpp $(COMMON_DIR)/wigner_compat.cpp
INCLUDES =

ifneq ($(strip $(EIGEN_DIR)),)
INCLUDES += -I$(EIGEN_DIR)
endif

all: main

main: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCES) -o $@

clean:
	rm -f main

# Legacy Unix helper. Prefer the repository root CMake build for cross-platform use.
