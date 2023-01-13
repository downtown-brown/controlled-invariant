CC=g++
CFLAGS=-O3

OBJECTS = \
	main.o \
	polyhedra.o \

%.o: %.cc
	$(CC) -c -o $@ $< $(CFLAGS)

debug: $(OBJECTS)
	$(CC) -o polygons $(OBJECTS) -lppl -lgmp

polygons: $(OBJECTS)
	$(CC) -o polygons $(OBJECTS) -lppl -lgmp
