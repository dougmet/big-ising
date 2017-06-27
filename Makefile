
COMMIT := $(shell git rev-parse HEAD)
BRANCH := $(shell git rev-parse --symbolic-full-name --abbrev-ref HEAD)
CC=g++
CFLAGS=-Wall -DPNG_DUMP -DCOMMIT=\"$(COMMIT)\" -DBRANCH=\"$(BRANCH)\"
LDFLAGS=
SOURCES=lodepng.cpp

all: ising isingpng isingmov

git:
	echo commit is $(COMMIT), branch is $(BRANCH)

ising: $(SOURCES) ising64.cpp
	$(CC) $(SOURCES) ising64.cpp $(CFLAGS) $(LDFLAGS) -o ising

isingpng: $(SOURCES) ising-png.cpp
	$(CC) $(SOURCES) ising-png.cpp $(CFLAGS) $(LDFLAGS) -o isingpng

isingmov: $(SOURCES) mov_huge_ising.cpp
	$(CC) $(SOURCES) mov_huge_ising.cpp $(CFLAGS) $(LDFLAGS) -o isingmov

clean:
	rm ising isingpng isingmov
