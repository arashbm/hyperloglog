CC = g++
CXX = g++
CXXFLAGS = -Werror -Wall -Wextra -O2 -funroll-loops -ffast-math -ftree-vectorize -mtune=native -std=c++11 -Wconversion -Wno-c++14-extensions -g

MurmurHash3.o: CXXFLAGS+=-Wno-implicit-fallthrough -Wno-sign-conversion

hyperloglog.hpp: hyperloglog.tpp
test_hll: MurmurHash3.o | hyperloglog.hpp
record_biases: MurmurHash3.o | hyperloglog.hpp 

clean:
	rm -f MurmurHash3.o test_hll record_biases
