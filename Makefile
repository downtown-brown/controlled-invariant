CC=g++
CFLAGS=-O3 -Wall -Werror -g -std=c++20
SRCDIR=src
OBJDIR=obj
DATADIR=data

EXE=invariant

SRC=$(wildcard $(SRCDIR)/*.cc)
OBJ=$(SRC:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) | $(DATADIR)
	$(CC) $(CFLAGS) $^ -lppl -lgmp -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR) $(DATADIR):
	mkdir -p $@

clean:
	@rm -rv $(OBJDIR)
