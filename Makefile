
CC=g++
CFLAGS=-c -Wall -DPNG_DUMP -I/opt/X11/include -I/opt/X11/include/freetype2 -I/Users/douglas/Cprogs/build/include
LDFLAGS=-llibpng
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
