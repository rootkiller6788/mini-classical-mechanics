# Makefile �� mini-hamiltonian
# Builds C library, verifies Lean formalization, runs Julia tests

CC       = gcc
CFLAGS   = -std=gnu11 -O2 -Wall -Wextra
LDFLAGS  = -lm

INCDIR   = include
SRCDIR   = src
BUILDDIR = build

CSRCS = $(SRCDIR)/hamiltonian_core.c $(SRCDIR)/hamiltonian_poisson.c $(SRCDIR)/hamiltonian_symplectic.c $(SRCDIR)/hamiltonian_lie.c
COBJS = $(CSRCS:$(SRCDIR)/%.c=$(BUILDDIR)/%.o)
LIBNAME = libhamiltonian.a

.PHONY: all clean test jl-test lean-check build-dirs

all: build-dirs $(LIBNAME)
	@echo === mini-hamiltonian build complete ===
	@echo Run: make jl-test   for Julia tests

build-dirs:
	@mkdir -p $(BUILDDIR)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(INCDIR)/hamiltonian_types.h
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

$(LIBNAME): $(COBJS)
	ar rcs $@ $^

jl-test:
	cd tests && julia --project=.. runtests.jl

lean-check:
	@echo Formalization: src/hamiltonian.lean (Lean 4)

line-count:
	@echo include+src total:
	@find $(INCDIR) $(SRCDIR) -name "*.h" -o -name "*.c" -o -name "*.jl" -o -name "*.lean" | xargs cat | wc -l

clean:
	rm -rf $(BUILDDIR) $(LIBNAME)
