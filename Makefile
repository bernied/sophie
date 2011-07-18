CFLAGS=-g -O2 -Wall -DNDEBUG -std=c99 $(OPTFLAGS)
OPTFLAGS=
LIBS=
SOURCES=sophie.c parse_cl.c
OBJECTS=sophie.o parse_cl.o

all: sophie

sophie: $(OBJECTS)
	$(CC) $(CFLAGS) $(LIBS) $(SOURCES) -o sophie

clean:
	rm -f *.o sophie
