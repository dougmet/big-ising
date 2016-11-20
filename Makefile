
CC=g++
CFLAGS=-Wall -DPNG_DUMP
LDFLAGS=
SOURCES=lodepng.cpp ising64.cpp

all: ising

ising:
	$(CC) $(SOURCES) $(CFLAGS) $(LDFLAGS) -o ising

