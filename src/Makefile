CC = g++
CFLAGS = -Wall -Werror -Wextra -g
STD = -std=c++20 
GCOV_FLAGS = -fprofile-arcs -ftest-coverage
GTEST_FLAGS = -lgtest
CLANG_FORMAT = clang-format-18
SRC_FILES = s21_matrix_oop.cc matrix_helpers_functions.cc
SRC_TEST = ./tests/tests.cc
OBJ_FILES = ./*.o
LOGS_DIR = ./logs
IGNORE = '/usr/include/*'

ifeq ($(shell uname), Linux)
	GTEST_FLAGS += -lgtest
	LCOVFLAG = 
else
	GTEST_FLAGS += -ld64
	CLANG_FORMAT = clang-format
	LCOVFLAG = --ignore-errors inconsistent,inconsistent
	IGNORE = '/usr/local/include/*'
endif

all: clean test s21_matrix_oop.a

s21_matrix_oop.a: 
	$(CC) $(CFLAGS) -c $(SRC_FILES)
	ar rcs s21_matrix_oop.a $(OBJ_FILES)
	rm -rf $(OBJ_FILES)

test: s21_matrix_oop.a $(SRC_TEST)
	$(CC) $(CFLAGS) $(STD) $(SRC_TEST) -L. -o s21_test s21_matrix_oop.a $(GTEST_FLAGS)
	./s21_test

gcov_report: 
	@echo "=== CLEANING ==="
	@rm -rf $(OBJ_FILES) s21_matrix_oop.a s21_test *.out $(LOGS_DIR) *.gcda *.gcno gcov_test
	
	@echo "=== COMPILING WITH COVERAGE ==="
	$(CC) $(CFLAGS) $(STD) --coverage $(SRC_FILES) $(SRC_TEST) -o gcov_test $(GTEST_FLAGS)
	
	@echo "=== RUNNING TESTS ==="
	./gcov_test
	
	@echo "=== GENERATING LCOV REPORT ==="
	@mkdir -p $(LOGS_DIR)
	@lcov $(LCOVFLAG) --capture --directory . --output-file $(LOGS_DIR)/coverage.info
	@lcov $(LCOVFLAG) --remove $(LOGS_DIR)/coverage.info $(IGNORE) --output-file $(LOGS_DIR)/coverage_filtered.info
	@genhtml $(LOGS_DIR)/coverage_filtered.info --output-directory $(LOGS_DIR)
	@open $(LOGS_DIR)/index.html

gcovr_report:
	@if ! command -v gcovr > /dev/null 2>&1; then \
		echo "Installing gcovr..."; \
		pip install gcovr; \
	fi
	
	@echo "=== CLEANING ==="
	@rm -rf $(OBJ_FILES) s21_matrix_oop.a s21_test *.out $(LOGS_DIR) *.gcda *.gcno gcov_test
	
	@echo "=== COMPILING WITH COVERAGE ==="
	$(CC) $(CFLAGS) $(STD) --coverage $(SRC_FILES) $(SRC_TEST) -o gcov_test $(GTEST_FLAGS)
	
	@echo "=== RUNNING TESTS ==="
	./gcov_test
	
	@echo "=== GENERATING GCOVR REPORT ==="
	@gcovr --root . --html --html-details -o $(LOGS_DIR)/coverage_report.html
	@open $(LOGS_DIR)/coverage_report.html

gcov_simple:
	@echo "=== CLEANING ==="
	@rm -rf $(OBJ_FILES) s21_matrix_oop.a s21_test *.out $(LOGS_DIR) *.gcda *.gcno gcov_test
	
	@echo "=== COMPILING WITH COVERAGE ==="
	$(CC) $(CFLAGS) $(STD) --coverage $(SRC_FILES) $(SRC_TEST) -o gcov_test $(GTEST_FLAGS)
	
	@echo "=== RUNNING TESTS ==="
	./gcov_test
	
	@echo "=== GENERATING COVERAGE SUMMARY ==="
	@gcov *.gcda
	
	@echo ""
	@echo "=== COVERAGE RESULTS ==="
	@for file in $(SRC_FILES); do \
		if [ -f "$$file.gcov" ]; then \
			echo "File: $$file"; \
			grep -E "^(Lines executed:|No lines)" "$$file.gcov" 2>/dev/null || echo "  Coverage data available"; \
			echo "---"; \
		fi; \
	done

valgrind: 
	valgrind --tool=memcheck --leak-check=yes ./tests/tests.cc

clean:
	rm -rf $(OBJ_FILES) 
	rm -rf s21_matrix_oop.a
	rm -rf s21_test
	rm -rf *.out
	rm -rf $(LOGS_DIR)
	rm -rf *.gcda *.gcno gcov_test
	# rm -rf rm ./tests/tests_log.log

rebuild: clean s21_matrix_oop.a test 

style:
	cp ../materials/linters/.clang-format ./
	cp ../materials/linters/.clang-format ./tests/
	$(CLANG_FORMAT) -n ./tests/*.cc
	$(CLANG_FORMAT) -i ./tests/*.cc
	$(CLANG_FORMAT) -n ./tests/*.cc
	$(CLANG_FORMAT) -n ./*.cc
	$(CLANG_FORMAT) -i ./*.cc
	$(CLANG_FORMAT) -n ./*.cc
	$(CLANG_FORMAT) -n ./*.h
	$(CLANG_FORMAT) -i ./*.h
	$(CLANG_FORMAT) -n ./*.h
	rm .clang-format
	rm ./tests/.clang-format