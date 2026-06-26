# Makefile — mini-celestial-mechanics
# Builds celestial mechanics C library and test suite.
#
# Usage:
#   make          — build library and tests
#   make test     — build and run tests
#   make clean    — remove build artifacts

CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -O2 -g -Iinclude
LDFLAGS  = -lm

SRCDIR   = src
INCDIR   = include
TESTDIR  = tests
BUILDDIR = build

SRCS     = $(wildcard $(SRCDIR)/*.c)
OBJS     = $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SRCS))
TESTS    = $(TESTDIR)/test_celestial.c
TESTBIN  = $(BUILDDIR)/test_celestial

.PHONY: all lib test clean

all: lib test

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TESTBIN): $(OBJS) $(TESTS) | $(BUILDDIR)
	$(CC) $(CFLAGS) $(TESTS) $(OBJS) -o $(TESTBIN) $(LDFLAGS)

lib: $(OBJS)

test: $(TESTBIN)
	./$(TESTBIN)

clean:
	rm -rf $(BUILDDIR)
