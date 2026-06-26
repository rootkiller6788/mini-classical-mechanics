CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c11 -Iinclude
LDFLAGS = -lm
SRCDIR = src
INCDIR = include
TESTSDIR = tests
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:.c=.o)

.PHONY: all clean test check lines

all: libvariational.a

libvariational.a: $(OBJS)
	ar rcs $@ $^

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TESTSDIR)/test_all $(TESTSDIR)/test_all.c -L. -lvariational $(LDFLAGS)
	./$(TESTSDIR)/test_all

check: all
	@echo "All source files compiled successfully."
	@echo "Library libvariational.a built."

lines:
	@echo "Headers:" && wc -l $(INCDIR)/*.h
	@echo "Sources:" && wc -l $(SRCDIR)/*.c
	@echo "Total:" && cat $(INCDIR)/*.h $(SRCDIR)/*.c | wc -l

clean:
	rm -f $(SRCDIR)/*.o libvariational.a $(TESTSDIR)/test_all
