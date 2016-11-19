
CC=g++
CFLAGS=-c -Wall -DPNG_DUMP -I/anaconda/include -I/anaconda/include/freetype2 -I/Users/douglas/Cprogs/build/include
LDFLAGS=-L/anaconda/lib -L/Users/douglas/Cprogs/build/lib -llibpng -llibpngwriter
SOURCES=ising64.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXE=ising

all: $(SOURCES) $(EXE)
    
$(EXE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm *.o ising
