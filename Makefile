all: test

CXXFLAGS=-Wall -pedantic -std=c++17 -O3 -I/usr/include/eigen3

test: test.o curves.o optimizer.o
	g++ -o $@ $^
