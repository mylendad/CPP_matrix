#include <cstring>
#include <iostream>
#include <stdexcept>

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_ = nullptr;

 public:
  S21Matrix(int rows, int cols);

  S21Matrix(int rows);

  S21Matrix();

  ~S21Matrix();

  S21Matrix(const S21Matrix& other);

  S21Matrix(S21Matrix&& other);

  void SumMatrix(const S21Matrix& other);

  S21Matrix& operator+(const S21Matrix& other);

  S21Matrix& operator+=(const S21Matrix& other);

  S21Matrix& operator-=(const S21Matrix& other);

  double& operator()(int i, int j);

  // double&& operator()(int i, int j);

  void copy(const S21Matrix& other);

  void cleanMatrix();

  const S21Matrix& operator=(const S21Matrix& left);

  void SubMatrix(const S21Matrix& other);

  int EqMatrix(const S21Matrix& other);

  int getRows() const;

  int getCols() const;

  void setRows(int rows);

  void setCols(int cols);

  double getMatrix(int row, int col) const;

  void setNumToMatrix(int row, int col, double number);

  void printMatr();
};

double s21_trimer_numb(double src);

S21Matrix::S21Matrix(int rows, int cols) : rows_{rows}, cols_{cols} {
  if (this->rows_ != 0 && this->cols_ != 0) {
    this->matrix_ = new double*[rows];
    for (int i = 0; i < rows; i++) {
      this->matrix_[i] = new double[cols]{0};
    }
  } else {
    throw std::exception();
  }
}

S21Matrix::S21Matrix(int rows) : S21Matrix(rows, rows) {}

S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

S21Matrix::~S21Matrix() { this->cleanMatrix(); }

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  if (this->rows_ != 0 && this->cols_ != 0) {
    this->copy(other);
  } else {
    throw std::exception();
  }
  // std::cerr << this->matrix_ << "!!!in constructor" << std::endl;
  // std::cerr << other.matrix_ << "!!!in constructor" << std::endl;
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_) {
  this->matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
  std::cerr << other.matrix_ << "!!!in perenos" << std::endl;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_ ||
      this->matrix_ == nullptr || other.matrix_ == nullptr) {
    throw std::exception();
  } else {
    S21Matrix(this->rows_, this->cols_);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        this->matrix_[i][j] += other.matrix_[i][j];
      }
    }
  }
}

const S21Matrix& S21Matrix::operator=(const S21Matrix& left) {
  if (&left == this) return *this;
  if (this->matrix_ != nullptr) {
    this->cleanMatrix();
  }
  this->copy(left);
  // std::cerr << this->matrix_ << "!!!inoperator=" << std::endl;
  // std::cerr << left.matrix_ << "!!!inoperator=" << std::endl;
  return *this;
}

S21Matrix& S21Matrix::operator+(const S21Matrix& other) {
  // creating result matrix

  this->SumMatrix(other);
  // std::cerr << this->matrix_ << "!!!inoperator+" << std::endl;
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  if (this->matrix_ != nullptr) {
    this->cleanMatrix();
  }
  this->copy(other);
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  if (this->matrix_ != nullptr) {
    this->cleanMatrix();
  }
  this->copy(other);
  this->SubMatrix(other);
  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (i < this->rows_ && i >= 0 && j < this->cols_ && j >= 0)
    return this->matrix_[i][j];
  else
    throw std::exception();
}

// double&& S21Matrix::operator()(int i, int j) {
//   this->matrix_[i][j];
//   this->matrix_ = nullptr;
//   this->rows_ = 0;
//   this->cols_ = 0;
//  }

void S21Matrix::copy(const S21Matrix& other) {
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  // std::cerr << this->matrix_ << "!!!copy begin" << std::endl;
  if (other.rows_ != 0 && other.cols_ != 0) {
    this->matrix_ = new double* [other.rows_] { 0 };
    for (int i = 0; i < other.rows_; i++) {
      // std::cerr << this->cols_ << "!!!this->cols_" << std::endl;
      // std::cerr << other.cols_ << "!!other.cols_" << std::endl;
      // std::cerr << "!!!new begin" << std::endl;
      this->matrix_[i] = new double[other.cols_]{0};
      // std::cerr << this->matrix_ << "!!!new end" << std::endl;
      std::memcpy(this->matrix_[i], other.matrix_[i],
                  other.cols_ * sizeof(double));
    }
  } else
    throw std::exception();

  // std::cerr << this->matrix_ << "!!!copy end" << std::endl;
  // return *this;
}

void S21Matrix::cleanMatrix() {
  for (int i = 0; i < this->rows_; i++) {
    delete[] this->matrix_[i];
  }
  delete[] this->matrix_;
  this->matrix_ = nullptr;
  this->rows_ = 0;
  this->cols_ = 0;
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::exception();
  } else if (this->matrix_ == nullptr || other.matrix_ == nullptr) {
    throw std::exception();
  } else {
    S21Matrix(this->rows_, this->cols_);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        this->matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  }
}

int S21Matrix::EqMatrix(const S21Matrix& other) {
  int result = 0;
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::exception();

  } else if (this->matrix_ == nullptr || other.matrix_ == nullptr) {
    throw std::exception();
  }
  S21Matrix(this->rows_, this->cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_ && result != 1; j++) {
      if (s21_trimer_numb(this->matrix_[i][j]) ==
          s21_trimer_numb(other.matrix_[i][j])) {
        result = 0;
        // std::cerr << "!!!result = 0" << std::endl;
      } else {
        // std::cerr << "!!!result = 1" << std::endl;
        // std::cerr << i << " " << j << std::endl;
        result = 1;
        // break;  // delete
      }
    }
  }

  return result;
}

int S21Matrix::getRows() const { return this->rows_; }

int S21Matrix::getCols() const { return this->cols_; }

void S21Matrix::setRows(int rows) {
  this->rows_ = rows;
  S21Matrix(this->rows_, this->cols_);
}

void S21Matrix::setCols(int cols) {
  this->rows_ = cols;
  S21Matrix(this->rows_, this->cols_);
}

double S21Matrix::getMatrix(int row, int col) const {
  return this->matrix_[row][col];
}

void S21Matrix::setNumToMatrix(int row, int col, double number) {
  this->matrix_[row][col] = number;
}

void S21Matrix::printMatr() {
  std::cerr << "!!!address:" << this->matrix_ << "!!!" << std::endl;
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

double s21_trimer_numb(double src) {
  double result = 0;
  char buffer[20];  // add define

  sprintf(buffer, "%.7lf", src);
  sscanf(buffer, "%lf", &result);
  return result;
}

int main() {
  try

  {
    S21Matrix pM(5, 4);

    S21Matrix pM2(5, 4);  //
    // S21Matrix* pMres = new S21Matrix(4, 5);
    // pM->SumMatrix->matrix_[0][0] = 1.0;
    pM2.setNumToMatrix(0, 0, 1.0);
    pM.SumMatrix(pM2);
    std::cout << pM.getMatrix(0, 0) << std::endl;
    S21Matrix pMres{pM};  // copy
    pMres.printMatr();
    std::cout << pMres.getMatrix(0, 0) << std::endl;
    pM.SubMatrix(pM2);
    std::cout << pM.getMatrix(0, 0) << std::endl;
    S21Matrix pM21 = pM2;
    pM21.SumMatrix(pM2);
    std::cout << pM21.getMatrix(0, 0) << std::endl;
    pM21.printMatr();
    pM2.printMatr();
    std::cerr << "!!!pMres eq pM2" << std::endl;
    std::cerr << pMres.EqMatrix(pM2) << std::endl;
    S21Matrix matrix2(std::move(pM2));
    std::cerr << "!!!pM21 eq matrix2" << std::endl;
    std::cerr << pM21.EqMatrix(matrix2) << std::endl;
    matrix2.printMatr();
    pM2.printMatr();
    pM2.printMatr();
    pM2.printMatr();
    S21Matrix MSum(5, 4);
    MSum.printMatr();
    MSum = (matrix2 + matrix2) + matrix2;
    std::cerr << "!!!matrix2" << std::endl;
    matrix2.printMatr();
    std::cerr << "!!!MSum" << std::endl;
    MSum.printMatr();
    MSum += matrix2;
    std::cerr << "!!!MSum += matrix2;" << std::endl;
    MSum.printMatr();
    std::cerr << MSum(0, 0) << std::endl;
    MSum(0, 0) = 777;
    std::cerr << MSum(0, 0) << std::endl;
    std::cerr << "!!!pM" << std::endl;
    pM.printMatr();
    std::cerr << "!!! MSum.EqMatrix(matrix2)" << std::endl;
    std::cerr << MSum.EqMatrix(matrix2) << std::endl;
    std::cerr << "!!! MSum.EqMatrix(pM)" << std::endl;
    std::cerr << MSum.EqMatrix(pM) << std::endl;
  } catch (const std::exception& err) {
    std::cerr << "[exception] " << err.what() << std::endl;
  }

  return 0;
}