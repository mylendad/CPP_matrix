#include "s21_matrix_oop.h"

#include <math.h>

#include <cstring>
#include <iostream>
#include <stdexcept>

S21Matrix::S21Matrix(int rows, int cols)
    : rows_{rows}, cols_{cols}, matrix_{nullptr} {
  if (rows_ <= 0 || cols_ <= 0) {
    throw std::invalid_argument("Matrix dimensions must be positive.");
  }

  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    try {
      matrix_[i] = new double[cols_]{0};
    } catch (std::bad_alloc &ba) {
      for (int j = 0; j < i; j++) {
        delete[] matrix_[j];
      }
      delete[] matrix_;
      matrix_ = nullptr;
      throw;
    }
  }
}

S21Matrix::S21Matrix(int rows) : S21Matrix(rows, rows) {}

S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

S21Matrix::~S21Matrix() { this->CleanMatrix(); }

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  this->CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_) {
  this->matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::invalid_argument("Different matrix dimensions.");
  }
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::invalid_argument("Different matrix dimensions.");
  }
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

bool S21Matrix::EqMatrix(const S21Matrix &other) noexcept {
  bool result = true;
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_ && result != false; j++) {
      if (std::abs((this->matrix_[i][j]) - (other.matrix_[i][j])) >= EPSILON) {
        result = false;
      } else {
        result = true;
      }
    }
  }
  return result;
}

void S21Matrix::MulNumber(const double number) noexcept {
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] *= number;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (this->rows_ != other.cols_ || this->cols_ != other.rows_) {
    throw std::invalid_argument(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix.");
  }

  S21Matrix res(this->rows_, other.cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < other.rows_; k++)
        res.matrix_[i][j] += (this->matrix_[i][k] * other.matrix_[k][j]);
    }
  }
  *this = res;
}

S21Matrix S21Matrix::Transpose() noexcept {
  S21Matrix res(this->cols_, this->rows_);
  for (int i = 0; i < this->cols_; i++) {
    for (int j = 0; j < this->rows_; j++) {
      res.matrix_[i][j] = this->matrix_[j][i];
    }
  }
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (this->rows_ != this->cols_) {
    throw std::invalid_argument("The matrix is not square.");
  }
  S21Matrix res(this->rows_);
  if (this->rows_ == 1) {
    this->matrix_[0][0] = this->matrix_[0][0];
  } else if (this->rows_ == 2) {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->rows_; j++) {
        res.matrix_[this->rows_ - 1 - i][this->rows_ - 1 - j] =
            this->matrix_[i][j] * pow(-1, ((i + 1) + (j + 1)));
      }
    }
  } else if (this->rows_ != 0) {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->rows_; j++) {
        res.matrix_[i][j] =

            CalcMinors(*this, i, j) *
            pow(-1, ((this->rows_ + 1 - i) + (this->rows_ + 1 - j)));
      }
    }
  }
  return res;
}

double S21Matrix::Determinant() {
  if (this->rows_ != this->cols_) {
    throw std::invalid_argument("The matrix is not square.");
  }
  double result = 0;
  if (this->rows_ == 1)
    result = this->matrix_[0][0];
  else {
    if (this->rows_ == 2)
      result = EndUnit(*this);
    else {
      for (int i = 0; i < this->rows_; i++) {
        double minor = CalcMinors(*this, 0, i) * this->matrix_[0][i];
        if (i % 2 > 0) minor *= -1;
        result += minor;
      }
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determin = this->Determinant();
  if (std::abs(determin) < EPSILON)
    throw std::invalid_argument("Matrix determinant is 0.");
  S21Matrix res(this->rows_);
  if (this->rows_ == 1) {
    res.matrix_[0][0] = 1 / this->matrix_[0][0];
  } else {
    S21Matrix minor_matrix = this->CalcComplements();
    S21Matrix transpose_matrix = minor_matrix.Transpose();
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->rows_; j++)
        res.matrix_[i][j] = (1 / determin) * transpose_matrix.matrix_[i][j];
    }
  }
  return res;
}

const S21Matrix &S21Matrix::operator=(const S21Matrix &other) noexcept {
  if (&other == this) return *this;
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->CopyMatrixSize(other);
  this->CopyMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (this == &other) return *this;
  CleanMatrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;

  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
  return *this;
}

bool S21Matrix::operator==(S21Matrix &other) noexcept {
  bool result = false;
  this->EqMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double number) const {
  S21Matrix result(*this);
  result.MulNumber(number);
  return result;
}

S21Matrix operator*(double number, const S21Matrix &matrix) {
  S21Matrix result(matrix);
  result.MulNumber(number);
  return result;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double number) {
  this->MulNumber(number);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  this->MulMatrix(other);
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  if (0 > i || this->rows_ < i + 1 || this->cols_ < j + 1 || 0 > j)
    throw std::out_of_range("Index is outside the matrix.");
  return this->matrix_[i][j];
}