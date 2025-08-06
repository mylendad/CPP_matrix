CC = g++
CFLAGS = -Wall -Werror -Wextra -g
STD = -std=c++20 
CHECK = -lcheck
GCOV_FLAGS = -fprofile-arcs -ftest-coverage
CLANG_FORMAT = clang-format-18
SRC_FILES = matrix.cc matrix_helpers_functions.cc
# SRC_TEST = ./tests/tests.c
OBJ_FILES = ./*.o

# ifeq ($(shell uname), Linux)
# 	CHECK += -lm -lsubunit
# else
# 	CHECK = -lcheck
# 	CLANG_FORMAT = clang-format
# endif

all: all, clean, test, s21_matrix_oop.a

s21_matrix_oop.a: 
	$(CC) $(CFLAGS) -c $(SRC_FILES)
	ar rcs s21_matrix_oop.a $(OBJ_FILES)
	rm -rf $(OBJ_FILES)

test: s21_matrix_oop.a
	$(CC) $(CFLAGS) $(SRC_TEST) -L. $(SRC_FILES) $(CHECK) -o s21_test 
	./s21_test

valgrind: 
	valgrind --tool=memcheck --leak-check=yes ./s21_test

	
clean:
	rm -rf $(OBJ_FILES) 
	rm -rf s21_matrix_oop.a
	rm -rf s21_test
	rm -rf *.out
	# rm -rf rm ./tests/tests_log.log

rebuild: clean s21_matrix_oop.a test 

style:
	# cp ../materials/linters/.clang-format ./
	# cp ../materials/linters/.clang-format ./tests/
	# $(CLANG_FORMAT) -n ./tests/*.c
	# $(CLANG_FORMAT) -i ./tests/*.c
	# $(CLANG_FORMAT) -n ./tests/*.c
	# $(CLANG_FORMAT) -n ./*.c
	# $(CLANG_FORMAT) -i ./*.c
	# $(CLANG_FORMAT) -n ./*.c
	# $(CLANG_FORMAT) -n ./*.h
	# $(CLANG_FORMAT) -i ./*.h
	# $(CLANG_FORMAT) -n ./*.h
	# rm .clang-format
	# rm ./tests/.clang-format
