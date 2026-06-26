# Makefile — mini-lagrangian Julia 模块
# make test: 运行所有测试
# make examples: 运行所有示例
# make clean: 清理临时文件

JULIA ?= julia
SRC_DIR = src
INCLUDE_DIR = include
TEST_DIR = tests
EXAMPLES_DIR = examples

.PHONY: all test examples clean check verify

all: test

test:
	@echo "=== Running mini-lagrangian tests ==="
	$(JULIA) $(TEST_DIR)/runtests.jl

examples:
	@echo "=== Running all examples ==="
	@for f in $(EXAMPLES_DIR)/*.jl; do \
		echo "--- $$f ---"; \
		$(JULIA) "$$f" || exit 1; \
	done

safety:
	@echo "=== Safety check: filler scan ==="
	@! grep -rn '_fn[0-9]' $(SRC_DIR)/ $(INCLUDE_DIR)/ 2>/dev/null || (echo "FAIL: _fn pattern" && exit 1)
	@! grep -rn '_aux[0-9]' $(SRC_DIR)/ $(INCLUDE_DIR)/ 2>/dev/null || (echo "FAIL: _aux pattern" && exit 1)
	@! grep -rn '_ext[0-9]' $(SRC_DIR)/ $(INCLUDE_DIR)/ 2>/dev/null || (echo "FAIL: _ext pattern" && exit 1)
	@! grep -rn 'TODO\|FIXME\|placeholder\|stub' $(SRC_DIR)/ $(INCLUDE_DIR)/ $(TEST_DIR)/ $(EXAMPLES_DIR)/ 2>/dev/null || (echo "FAIL: stub markers" && exit 1)
	@echo "=== Safety check PASSED ==="

check: safety test

verify:
	@echo "=== Verifying COMPLETE status ==="
	@lines=$$(cat $(SRC_DIR)/*.jl $(INCLUDE_DIR)/*.jl 2>/dev/null | wc -l); \
	echo "include/ + src/ total lines: $$lines"; \
	if [ "$$lines" -ge 3000 ]; then \
		echo "PASS: Line count >= 3000"; \
	else \
		echo "FAIL: Line count $$lines < 3000"; exit 1; \
	fi

clean:
	@rm -f *.jl.cov *.mem
	@echo "Cleaned."
