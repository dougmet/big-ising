
COMMIT := $(shell git rev-parse HEAD)
BRANCH := $(shell git rev-parse --symbolic-full-name --abbrev-ref HEAD)
CC=g++
CFLAGS=-O3 -Wall -DCOMMIT=\"$(COMMIT)\" -DBRANCH=\"$(BRANCH)\"
LDFLAGS=
SOURCES=src/lodepng.cpp

# -DPNG_DUMP

all: ising isingpng isingmov

git:
	echo commit is $(COMMIT), branch is $(BRANCH)

ising: $(SOURCES) src/ising64.cpp
	$(CC) $(SOURCES) src/ising64.cpp $(CFLAGS) $(LDFLAGS) -o ising

isingpng: $(SOURCES) src/ising-png.cpp
	$(CC) $(SOURCES) src/ising-png.cpp $(CFLAGS) $(LDFLAGS) -o isingpng

isingmov: $(SOURCES) src/mov_huge_ising.cpp
	$(CC) $(SOURCES) src/mov_huge_ising.cpp $(CFLAGS) $(LDFLAGS) -o isingmov

isingclusters: $(SOURCES) src/ising-clusters.cpp
	$(CC) $(SOURCES) src/ising-clusters.cpp -DPNG_DUMP $(CFLAGS) $(LDFLAGS) -o isingclusters

clean:
	rm ising isingpng isingmov
