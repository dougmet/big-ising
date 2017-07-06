
COMMIT := $(shell git rev-parse HEAD)
BRANCH := $(shell git rev-parse --symbolic-full-name --abbrev-ref HEAD)
CC=g++
CFLAGS=-O3 -Wall -DCOMMIT=\"$(COMMIT)\" -DBRANCH=\"$(BRANCH)\"
LDFLAGS=
SOURCES=lodepng.cpp

# -DPNG_DUMP

all: ising isingpng isingmov

git:
	echo commit is $(COMMIT), branch is $(BRANCH)

ising: $(SOURCES) ising64.cpp
	$(CC) $(SOURCES) ising64.cpp $(CFLAGS) $(LDFLAGS) -o ising

isingpng: $(SOURCES) ising-png.cpp
	$(CC) $(SOURCES) ising-png.cpp $(CFLAGS) $(LDFLAGS) -o isingpng

isingmov: $(SOURCES) mov_huge_ising.cpp
	$(CC) $(SOURCES) mov_huge_ising.cpp $(CFLAGS) $(LDFLAGS) -o isingmov

isingclusters: $(SOURCES) ising-clusters.cpp
	$(CC) $(SOURCES) ising-clusters.cpp -DPNG_DUMP $(CFLAGS) $(LDFLAGS) -o isingclusters

clean:
	rm ising isingpng isingmov
