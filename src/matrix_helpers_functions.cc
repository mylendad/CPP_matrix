#include "s21_matrix_oop.h"

void S21Matrix::CopyMatrixSize(const S21Matrix& other) {
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
}

void S21Matrix::CopyMatrix(const S21Matrix& other) {
  if (this->matrix_ != nullptr) {
    for (int i = 0; i < this->rows_; i++) {
      delete[] this->matrix_[i];
    }
    delete[] this->matrix_;
  }
  if (other.rows_ != 0 && other.cols_ != 0) {
    this->matrix_ = new double* [this->rows_] { 0 };
    for (int i = 0; i < this->rows_; i++) {
      this->matrix_[i] = new double[this->cols_]{0};
      std::memcpy(this->matrix_[i], other.matrix_[i],
                  this->cols_ * sizeof(double));
    }
  } else
    throw std::exception();
}

void S21Matrix::CopyNumbersMatrix(const S21Matrix& other) {
  int min_rows = 0;
  int min_cols = 0;
  if (other.rows_ < this->rows_)
    min_rows = other.rows_;
  else
    min_rows = this->rows_;
  if (other.cols_ < this->cols_)
    min_cols = other.cols_;
  else
    min_cols = this->cols_;
  if (other.rows_ != 0 && other.cols_ != 0) {
    for (int i = 0; i < min_rows; i++) {
      std::memcpy(this->matrix_[i], other.matrix_[i],
                  min_cols * sizeof(double));
    }
  } else
    throw std::exception();
}

void S21Matrix::CleanMatrix() {
  if (this->matrix_ == nullptr) return;
  for (int i = 0; i < this->rows_; i++) {
    if (this->matrix_[i] != nullptr) {
      delete[] this->matrix_[i];
      this->matrix_[i] = nullptr;
    }
  }
  delete[] this->matrix_;
  this->matrix_ = nullptr;
  this->rows_ = 0;
  this->cols_ = 0;
}

int S21Matrix::GetRows() const { return this->rows_; }

int S21Matrix::GetCols() const { return this->cols_; }

void S21Matrix::SetRows(int rows) {
  S21Matrix temp(rows, this->cols_);
  temp.CopyNumbersMatrix(*this);
  *this = temp;
  this->rows_ = rows;
}

void S21Matrix::SetCols(int cols) {
  S21Matrix temp(this->rows_, cols_);
  temp.CopyNumbersMatrix(*this);
  *this = temp;
  this->cols_ = cols;
}

double S21Matrix::GetMatrix(int row, int col) const {
  return this->matrix_[row][col];
}

void S21Matrix::SetNumToMatrix(int row, int col, double number) {
  this->matrix_[row][col] = number;
}

void S21Matrix::PrintMatr() {
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      if (i == this->rows_ - 1 && j == this->cols_ - 1)
        std::cout << this->matrix_[i][j] << std::endl;
      else if (j == this->cols_ - 1)
        std::cout << this->matrix_[i][j] << std::endl;
      else
        std::cout << this->matrix_[i][j] << " ";
    }
  }
}

S21Matrix S21Matrix::CreateDeterminateMatrix(S21Matrix& other, int row,
                                             int col) {
  S21Matrix res(other.rows_ - 1);
  for (int i = 0, x = 0; i < other.rows_; i++) {
    for (int j = 0, y = 0; j < other.rows_; j++) {
      if (j != col && i != row) {
        res.matrix_[x][y] = other.matrix_[i][j];
        y++;
      }
      if (y == other.rows_ - 1) {
        y = 0;
        x++;
      }
    }
  }
  return res;
}

double S21Matrix::CalcMinors(S21Matrix& other, int row, int col) {
  double minor = 0.0;
  S21Matrix temp;
  if (other.rows_ > 2) temp = CreateDeterminateMatrix(other, row, col);
  if (temp.rows_ == 2) minor = EndUnit(temp);
  if (temp.rows_ > 2)
    for (int j = 0; j < temp.rows_ - 1; j++) {
      double a = CalcMinors(temp, 0, j) * temp.matrix_[0][j];
      if (!IsEvenNumber(j)) a *= -1;
      minor += a;
    }
  return minor;
}

double S21Matrix::EndUnit(S21Matrix& other) {
  double result = 0.0;
  result = (other.matrix_[0][0] * other.matrix_[1][1]) -
           (other.matrix_[0][1] * other.matrix_[1][0]);
  return result;
}

double S21Matrix::TrimerNumb(double src) {
  double result = 0;
  char buffer[BUFFER_TRIMMER];
  sprintf(buffer, "%.7lf", src);
  sscanf(buffer, "%lf", &result);
  return result;
}

bool S21Matrix::IsEvenNumber(int number) const {
  bool result = FALSE;
  if (number % 2 > 0) result = TRUE;
  return result;
}