
CC=g++
CFLAGS=-c -Wall
LDFLAGS=
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
