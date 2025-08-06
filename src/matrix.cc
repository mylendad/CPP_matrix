#include "matrix.hh"

#include <math.h>

#include <cstring>
#include <iostream>
#include <stdexcept>

S21Matrix::S21Matrix(int rows, int cols) : rows_{rows}, cols_{cols} {
  if (this->rows_ < 0 || this->cols_ < 0) {
    throw std::invalid_argument("The argument cannot have a negative value.");
  }
  try {
    this->matrix_ = new double*[rows];
  } catch (std::bad_alloc& ba) {
    std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
    return;
  }

  for (int i = 0; i < rows; i++) {
    try {
      this->matrix_[i] = new double[cols]{0};
    } catch (std::bad_alloc& ba) {
      delete[] this->matrix_;
      this->matrix_ = nullptr;
      std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
    }
  }
}

S21Matrix::S21Matrix(int rows) : S21Matrix(rows, rows) {}

S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

S21Matrix::~S21Matrix() { this->CleanMatrix(); }

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  this->CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_) {
  this->matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::invalid_argument("Different matrix dimensions.");
  }
  S21Matrix(this->rows_, this->cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::invalid_argument("Different matrix dimensions.");
  }
  S21Matrix(this->rows_, this->cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

int S21Matrix::EqMatrix(const S21Matrix& other) {
  int result = 0;
  S21Matrix(this->rows_, this->cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_ && result != 1; j++) {
      if (TrimerNumb(this->matrix_[i][j]) == TrimerNumb(other.matrix_[i][j])) {
        result = 0;
      } else {
        result = 1;
      }
    }
  }
  return result;
}

void S21Matrix::MulNumber(const double number) {
  S21Matrix(this->rows_, this->cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] *= number;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->rows_ != other.cols_ || this->cols_ != other.rows_) {
    throw std::invalid_argument(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix.");
  }

  S21Matrix(this->rows_, this->cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      for (int k = 0; k < other.rows_; k++)
        this->matrix_[i][j] += (this->matrix_[i][k] * other.matrix_[k][j]);
    }
  }
}

S21Matrix S21Matrix::Transpose() {
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
        if (!IsEvenNumber(i)) minor *= -1;
        result += minor;
      }
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determin = this->Determinant();
  if (determin == 0) throw std::invalid_argument("Matrix determinant is 0.");
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

const S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (&other == this) return *this;
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->CopyMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  if (&other == this) return *this;
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
  return *this;
}

S21Matrix& S21Matrix::operator+(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->CopyMatrix(other);
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->CopyMatrix(other);
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->CopyMatrix(other);
  this->MulMatrix(other);
  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (0 > i || this->rows_ < i + 1 || this->cols_ < i + 1 || 0 > j ||
      this->rows_ < j + 1 || this->cols_ < j + 1)
    throw std::out_of_range("Index is outside the matrix.");
  if (i < this->rows_ && i >= 0 && j < this->cols_ && j >= 0)
    return this->matrix_[i][j];
  else
    throw std::exception();
}

int main() {
  try

  {
    S21Matrix pM(2, 2);

    S21Matrix pM2(2, 2);  //
    pM2.SetNumToMatrix(0, 0, 1.0);
    pM2.PrintMatr();
    pM.PrintMatr();
    pM.SumMatrix(pM2);
    std::cout << pM.GetMatrix(0, 0) << std::endl;
    S21Matrix pMres{pM};  // CopyMatrix
    pMres.PrintMatr();
    std::cout << pMres.GetMatrix(0, 0) << std::endl;
    pM.SubMatrix(pM2);
    std::cout << pM.GetMatrix(0, 0) << std::endl;
    S21Matrix pM21 = pM2;
    pM21.SumMatrix(pM2);
    std::cout << pM21.GetMatrix(0, 0) << std::endl;
    pM21.PrintMatr();
    pM2.PrintMatr();
    std::cout << "!!!pMres eq pM2" << std::endl;
    std::cout << pMres.EqMatrix(pM2) << std::endl;
    S21Matrix matrix2(std::move(pM2));
    std::cout << "!!!pM21 eq matrix2" << std::endl;
    std::cout << pM21.EqMatrix(matrix2) << std::endl;
    matrix2.PrintMatr();
    pM2.PrintMatr();
    pM2.PrintMatr();
    pM2.PrintMatr();
    S21Matrix MSum(2, 2);
    MSum.PrintMatr();
    MSum = (matrix2 + matrix2) + matrix2;
    std::cout << "!!!matrix2" << std::endl;
    matrix2.PrintMatr();
    std::cout << "!!!MSum" << std::endl;
    MSum.PrintMatr();
    MSum += matrix2;
    std::cout << "!!!MSum += matrix2;" << std::endl;
    MSum.PrintMatr();
    std::cout << MSum(0, 0) << std::endl;
    MSum(0, 0) = 4.0;
    MSum(0, 1) = 3.0;
    // MSum(0, 2) = 3;
    // MSum(0, 3) = 4;
    MSum(1, 0) = 6.0;
    MSum(1, 1) = 3.0;
    // MSum(1, 2) = 2;
    // MSum(1, 3) = 8;
    // MSum(2, 0) = 5;
    // MSum(2, 1) = 2;
    // MSum(2, 2) = 1;
    // MSum(2, 3) = 12;
    // MSum(3, 0) = 13;
    // MSum(3, 1) = 14;
    // MSum(3, 2) = 15;
    // MSum(3, 3) = 16;
    std::cout << "!!!MSum" << std::endl;
    MSum.PrintMatr();
    S21Matrix pM22;
    pM22 = MSum.Transpose();  ////
    std::cout << "!!! MSum.Transpose();" << std::endl;
    std::cout << "!!!pM22" << std::endl;
    pM22.PrintMatr();
    std::cout << "!!! MSum.GetRows() " << MSum.GetRows() << std::endl;
    S21Matrix pM23 = MSum.CalcComplements();
    std::cerr << "!!!pM23 = MSum.CalcComplements();" << std::endl;
    std::cerr << "!!!pM23" << std::endl;
    pM23.PrintMatr();

    std::cerr << "!!!determ;" << std::endl;
    std::cerr << MSum.Determinant() << std::endl;
    // pM23.PrintMatr();

    MSum(0, 0) = 4.0;
    MSum(0, 1) = 7.0;
    // MSum(0, 2) = 3;
    // MSum(0, 3) = 4;
    MSum(1, 0) = 2.0;
    MSum(1, 1) = 6.0;
    MSum.PrintMatr();
    S21Matrix pM24 = MSum.InverseMatrix();
    std::cout << "!!! MSum.InverseMatrix();" << std::endl;
    pM24.PrintMatr();

    std::cout << MSum(0, 0) << std::endl;
    std::cout << MSum(0, 0) << std::endl;
    std::cout << "!!!pM" << std::endl;
    pM.PrintMatr();
    std::cout << "!!! MSum.EqMatrix(matrix2)" << std::endl;
    std::cout << MSum.EqMatrix(matrix2) << std::endl;
    std::cout << "!!! MSum.EqMatrix(pM)" << std::endl;
    std::cout << MSum.EqMatrix(pM) << std::endl;
  } catch (const std::exception& err) {
    std::cerr << "[exception] " << err.what() << std::endl;
  }

  return 0;
}