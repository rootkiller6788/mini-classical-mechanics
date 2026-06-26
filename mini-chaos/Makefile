# Makefile - mini-chaos C library
# Build: make, Test: make test, Clean: make clean

CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2 -I include
LDFLAGS  = -lm
AR       = ar
ARFLAGS  = rcs

SRCS     = src/chaos_types.c src/chaos_flows.c src/chaos_maps.c \
           src/chaos_lyapunov.c src/chaos_bifurcation.c src/chaos_fractals.c \
           src/chaos_embedding.c src/chaos_analysis.c
OBJS     = $(SRCS:.c=.o)
LIB      = libchaos.a

TEST_SRC = tests_c/test_chaos.c
TEST_BIN = tests_c/test_chaos

EXAMPLE_SRCS = examples_c/example_lorenz.c examples_c/example_logistic.c examples_c/example_fractals.c
EXAMPLE_BINS = $(EXAMPLE_SRCS:.c=)

.PHONY: all clean test examples

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIB) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) -L. -lchaos $(LDFLAGS) -o $(TEST_BIN)
	./$(TEST_BIN)

examples: $(LIB) $(EXAMPLE_BINS)

examples_c/example_lorenz: examples_c/example_lorenz.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lchaos $(LDFLAGS) -o $@

examples_c/example_logistic: examples_c/example_logistic.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lchaos $(LDFLAGS) -o $@

examples_c/example_fractals: examples_c/example_fractals.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lchaos $(LDFLAGS) -o $@

clean:
	rm -f $(OBJS) $(LIB) $(TEST_BIN) $(EXAMPLE_BINS)
	rm -f tests_c/*.o examples_c/*.o
