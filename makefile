CC = g++
CXX = g++
CXXFLAGS = -Werror -Wall -Wextra -O3 -funroll-loops -ffast-math -ftree-vectorize -mtune=native -std=c++11 -Wconversion -pthread
SHELL=bash

hyperloglog.hpp: hyperloglog.tpp
test_hll: hyperloglog.hpp MurmurHash3.o
record_biases: hyperloglog.hpp MurmurHash3.o

clean:
	rm -f MurmurHash3.o test_hll record_biases
