#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <math.h>

#include <cstring>
#include <iostream>
#include <stdexcept>

#define TRUE 1
#define FALSE 0
#define EPSILON 1e-7

class S21Matrix {
 public:
  S21Matrix(int rows, int cols);

  S21Matrix(int rows);

  S21Matrix();

  ~S21Matrix();

  S21Matrix(const S21Matrix &other);

  S21Matrix(S21Matrix &&other);

  void SumMatrix(const S21Matrix &other);

  void SubMatrix(const S21Matrix &other);

  bool EqMatrix(const S21Matrix &other) noexcept;

  void MulNumber(const double number) noexcept;

  void MulMatrix(const S21Matrix &other);

  S21Matrix Transpose() noexcept;

  S21Matrix CalcComplements();

  double Determinant();

  S21Matrix InverseMatrix();

  const S21Matrix &operator=(const S21Matrix &other) noexcept;

  S21Matrix &operator=(S21Matrix &&other) noexcept;

  bool operator==(S21Matrix &other) noexcept;

  S21Matrix operator+(const S21Matrix &other) const;

  S21Matrix operator-(const S21Matrix &other) const;

  S21Matrix operator*(const S21Matrix &other) const;

  S21Matrix operator*(const double number) const;

  friend S21Matrix operator*(double number, const S21Matrix &matrix);

  S21Matrix &operator+=(const S21Matrix &other);

  S21Matrix &operator-=(const S21Matrix &other);

  S21Matrix &operator*=(const S21Matrix &other);

  S21Matrix &operator*=(const double number);

  double &operator()(int i, int j);

  int GetRows() const;

  int GetCols() const;

  void SetRows(int rows);

  void SetCols(int cols);

  double GetMatrix(int row, int col) const;

  void SetNumToMatrix(int row, int col, double number) noexcept;

  void PrintMatr();

 private:
  int rows_, cols_;

  double **matrix_ = nullptr;

  void CopyMatrixSize(const S21Matrix &other) noexcept;

  void CopyMatrix(const S21Matrix &other);

  void CopyNumbersMatrix(const S21Matrix &other);

  void CleanMatrix() noexcept;

  double CalcMinors(S21Matrix &other, int row, int col) noexcept;

  S21Matrix CreateDeterminateMatrix(S21Matrix &other, int row,
                                    int col) noexcept;

  double EndUnit(S21Matrix &other) noexcept;

  bool IsEvenNumber(int number) const;
};

#endif  // S21_MATRIX_OOP_H