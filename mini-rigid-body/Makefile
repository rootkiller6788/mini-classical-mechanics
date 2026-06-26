# Makefile — Mini Rigid Body Mechanics
# Reference: Goldstein Ch.4-5, Landau Ch.6, MIT 8.012

CC      := gcc
CFLAGS  := -Wall -Wextra -std=c11 -O2 -g
LDFLAGS := -lm

INCDIR  := include
SRCDIR  := src
TESTDIR := tests

SRCS    := $(wildcard $(SRCDIR)/*.c)
OBJS    := $(SRCS:.c=.o)
TESTSRC := $(TESTDIR)/test_rigid_body.c
TESTEXE := $(TESTDIR)/test_rigid_body

.PHONY: all clean test

all: lib test

# Build object files from all src/*.c (excluding .lean files)
%.o: %.c
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

# Build and run tests
test: $(TESTEXE)
	./$(TESTEXE)

$(TESTEXE): $(OBJS) $(TESTSRC)
	$(CC) $(CFLAGS) -I$(INCDIR) $^ $(LDFLAGS) -o $@

clean:
	rm -f $(SRCDIR)/*.o $(TESTEXE) *.o
