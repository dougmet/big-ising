
CC=g++
CFLAGS=-Wall -DPNG_DUMP
LDFLAGS=
SOURCES=lodepng.cpp

all: ising isingpng isingmov

ising:
	$(CC) $(SOURCES) ising64.cpp $(CFLAGS) $(LDFLAGS) -o ising

isingpng:
	$(CC) $(SOURCES) ising-png.cpp $(CFLAGS) $(LDFLAGS) -o isingpng

isingmov:
	$(CC) $(SOURCES) mov_huge_ising.cpp $(CFLAGS) $(LDFLAGS) -o isingmov

clean:
	rm ising isingpng isingmov
