CC=g++
CFLAGS=-O3

SRCDIR=src
OBJDIR=obj
BINDIR=.

EXE=$(BINDIR)/polygons

SRC=$(wildcard $(SRCDIR)/*.cc)
OBJ=$(SRC:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJ)
	echo $(OBJ)
	$(CC) $^ -lppl -lgmp -o $@
