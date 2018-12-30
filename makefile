CC = g++
CXX = g++
CXXFLAGS = -Werror -Wall -Wextra -O2 -funroll-loops -ffast-math -ftree-vectorize -mtune=native -std=c++14 -Wconversion -g

all: MurmurHash3.o hll_tests

MurmurHash3.o: CXXFLAGS+=-Wno-implicit-fallthrough -Wno-sign-conversion

hyperloglog.hpp: hyperloglog.tpp
record_biases: MurmurHash3.o | hyperloglog.hpp
estimate_distribution: MurmurHash3.o | hyperloglog.hpp

hll_tests: MurmurHash3.o tests.o | hyperloglog.hpp

clean:
	rm -f MurmurHash3.o test_hll record_biases hll_tests.o
