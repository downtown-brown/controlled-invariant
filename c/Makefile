CC=g++
CFLAGS=-O3 -Ofast -Wall -Werror

SRCDIR=src
OBJDIR=obj
DATADIR=data

EXE=invariant

SRC=$(wildcard $(SRCDIR)/*.cc)
OBJ=$(SRC:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) | $(DATADIR)
	echo $(OBJ)
	$(CC) $^ -lppl -lgmp -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR) $(DATADIR):
	mkdir -p $@

clean:
	@rm -rv $(OBJDIR)
