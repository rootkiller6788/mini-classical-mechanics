# Makefile — mini-continuum (Julia + C)
# Julia tests and examples | C library for continuum types

CC = gcc
CFLAGS = -Wall -Wextra -O2 -Iinclude -std=c99
AR = ar

C_SRCS = $(wildcard src/*.c)
C_OBJS = $(C_SRCS:.c=.o)
LIB = build/libcontinuum.a

.PHONY: all test examples clean

all: $(LIB)
	@echo "mini-continuum C library built: $(LIB)"

$(LIB): $(C_OBJS)
	@mkdir -p build
	$(AR) rcs $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	@echo "=== Julia Tests ==="
	@which julia >/dev/null 2>&1 && cd tests && julia --project=.. runtests.jl || echo "Julia not found, skipping Julia tests"
	@echo "=== C Tests ==="
	$(CC) $(CFLAGS) tests/test_c.c -Lbuild -lcontinuum -lm -o build/test_c
	./build/test_c

examples: all
	@echo "=== Julia Examples ==="
	@which julia >/dev/null 2>&1 && cd examples && julia --project=.. wave_and_beam.jl || echo "Julia not found"
	@which julia >/dev/null 2>&1 && cd examples && julia --project=.. buckling.jl || true
	@which julia >/dev/null 2>&1 && cd examples && julia --project=.. seismic_waves.jl || true

clean:
	rm -rf build src/*.o

