# Makefile — mini-newtonian 构建与测试
# Julia (主) + C (双语言实现层)

.PHONY: all test examples clean cpp-build cpp-clean c-build c-test

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -march=native
LDFLAGS  = -lm
C_SRCDIR = src
C_INCDIR = include
C_BUILDDIR = build/c
C_SRCS   = $(wildcard $(C_SRCDIR)/*.c)
C_OBJS   = $(patsubst $(C_SRCDIR)/%.c,$(C_BUILDDIR)/%.o,$(C_SRCS))
C_TARGET = $(C_BUILDDIR)/libnewtonian.a
C_TEST   = $(C_BUILDDIR)/test_newtonian

all: c-build

# ── Julia 测试 ──
test:
	julia tests/runtests.jl

# ── C 编译 ──
$(C_BUILDDIR)/%.o: $(C_SRCDIR)/%.c
	@mkdir -p $(C_BUILDDIR)
	$(CC) $(CFLAGS) -I$(C_INCDIR) -c $< -o $@

$(C_TARGET): $(C_OBJS)
	@mkdir -p $(C_BUILDDIR)
	ar rcs $@ $^

c-build: $(C_TARGET)
	@echo "[C] libnewtonian.a built ($(words $(C_OBJS)) objects)"

# ── C 测试可执行文件 ──
$(C_TEST): tests/test_newtonian.c $(C_TARGET)
	@mkdir -p $(C_BUILDDIR)
	$(CC) $(CFLAGS) -I$(C_INCDIR) $< -L$(C_BUILDDIR) -lnewtonian $(LDFLAGS) -o $@

c-test: $(C_TEST)
	@echo "=== C Tests ==="
	@$(C_TEST)

# ── Julia 示例 ──
examples:
	julia examples/projectile.jl
	julia examples/pendulum.jl
	julia examples/oscillator.jl
	julia examples/central_force.jl
	julia examples/nbody.jl

example-%:
	julia examples/$*.jl

# ── C++ 加速库 ──
cpp-build:
	cd cpp && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build .

cpp-clean:
	rm -rf cpp/build

# ── 清理 ──
clean:
	rm -rf build
	rm -rf .julia
	rm -rf *.ji

quick:
	julia -e 'include("tests/runtests.jl")'

info:
	@echo "=== mini-newtonian ==="
	@echo "C headers (include/):"
	@ls -la include/
	@echo ""
	@echo "C sources (src/):"
	@ls -la src/*.c 2>/dev/null || echo "(none)"
	@echo ""
	@echo "Julia sources (src/):"
	@ls -la src/*.jl 2>/dev/null || echo "(none)"
	@echo ""
	@echo "Line count (include/ + src/.c):"
	@find . \\( -path "*/include/*.h" -o -path "*/src/*.c" \\) -exec cat {} + 2>/dev/null | wc -l
