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

  const S21Matrix& operator+(const S21Matrix& other);

  const S21Matrix copy(const S21Matrix& other);

  void cleanMatrix();

  const S21Matrix& operator=(const S21Matrix& left);

  void SubMatrix(const S21Matrix& other);

  void EqMatrix(const S21Matrix& other);

  int getRows() const { return this->rows_; }

  int getCols() const { return this->cols_; }

  void setRows(int rows);

  void setCols(int cols);

  double getMatrix(int row, int col) const {
    // for (int i = 0; i < this->rows_; i++) {
    //   for (int j = 0; j < this->cols_; j++) {
    return this->matrix_[row][col];
    //   }
    // }
  }

  void setNumToMatrix(int row, int col, double number);

  void printMatr();
};

S21Matrix::S21Matrix(int rows, int cols) : rows_{rows}, cols_{cols} {
  this->matrix_ = new double*[rows];
  for (int i = 0; i < rows; i++) {
    this->matrix_[i] = new double[cols]{0};
  }
}

S21Matrix::S21Matrix(int rows) : S21Matrix(rows, rows) {}

S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

S21Matrix::~S21Matrix() { this->cleanMatrix(); }

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  this->matrix_ = new double*[other.rows_];  // вынести
  for (int i = 0; i < other.rows_; i++) {
    matrix_[i] = new double[this->cols_]{0};
    std::memcpy(this->matrix_[i], other.matrix_[i],
                other.cols_ * sizeof(double));
    // std::cerr << "!" << *this->matrix_[i] << "!" << " !" <<
    // *other.matrix_[i]
    // << "!" << std::endl;
  }
  std::cerr << this->matrix_ << "!!!in constructor" << std::endl;
  std::cerr << other.matrix_ << "!!!in constructor" << std::endl;
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_) {
  this->matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
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

const S21Matrix& S21Matrix::operator+(const S21Matrix& other) {
  // creating result matrix

  this->SumMatrix(other);

  return *this;
}

const S21Matrix S21Matrix::copy(const S21Matrix& other) {
  // вынести
  for (int i = 0; i < other.rows_; i++) {
    this->matrix_[i] = new double[this->cols_]{0};
    std::memcpy(this->matrix_[i], other.matrix_[i],
                other.cols_ * sizeof(double));
  }
  return *this;
}

void S21Matrix::cleanMatrix() {
  for (int i = 0; i < this->rows_; i++) {
    delete[] this->matrix_[i];
  }
  delete[] this->matrix_;
}

const S21Matrix& S21Matrix::operator=(const S21Matrix& left) {
  if (&left == this) return *this;
  if (this->matrix_ != nullptr) {
    // for (int i = 0; i < this->rows_; i++) {
    //   delete[] this->matrix_[i];
    // }
    // delete[] this->matrix_;
    this->cleanMatrix();  // вернуть
  }
  this->copy(left);  // не срабатывает мой конструктор
  std::cerr << this->matrix_ << "!!!inoperator=" << std::endl;
  std::cerr << left.matrix_ << "!!!inoperator=" << std::endl;
  // испраить или оставить так?
  return *this;  // исправить? чтобы создавался новый адрес в памти как в
                 // конструкторе копирования
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_ ||
      this->matrix_ == nullptr || other.matrix_ == nullptr) {
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

void S21Matrix::EqMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::exception();

  } else if (this->matrix_ == nullptr || other.matrix_ == nullptr) {
    throw std::exception();
  }
  S21Matrix(this->rows_, this->cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

// S21Matrix::S21Matrix(S21Matrix&& other) {
//   if (rows_ * cols_ == other.rows_ * other.cols_) {
//     std::memcpy(p_, other.p_, other.cols_ * other.rows_ * sizeof(double));
//   } else {
//     delete[] p_;
//     p_ = new double[other.rows_ * other.cols_]();
//     std::memcpy(p_, other.p_, other.cols_ * other.rows_ * sizeof(double));
//   }
//   rows_ = other.rows_;
//   cols_ = other.cols_;
//   delete other.p_;
//   other.rows_ = 0;
//   other.cols_ = 0;
// }

// int S21Matrix::getRows() const { return this->rows_; }

// int S21Matrix::getCols() const { return this->cols_; }

void S21Matrix::setRows(int rows) {
  this->rows_ = rows;
  S21Matrix(this->rows_, this->cols_);
}

void S21Matrix::setCols(int cols) {
  this->rows_ = cols;
  S21Matrix(this->rows_, this->cols_);
}
// double** setNumToMatrix(int rows, int cols, double number,
//                    const S21Matrix& other) {
//   for (int i = 0; i < rows; i++) {
//     for (int j = 0; j < cols; j++) {
//       other.matrix_[i][j] = number;
//     }
//   }
// }

// double S21Matrix::getMatrix(int row, int col) const {
//   // for (int i = 0; i < this->rows_; i++) {
//   //   for (int j = 0; j < this->cols_; j++) {
//   return this->matrix_[row][col];
//   //   }
//   // }
// }

void S21Matrix::setNumToMatrix(int row, int col, double number) {
  // for (int i = 0; i < this->rows_; i++) {
  //   for (int j = 0; j < this->cols_; j++) {
  this->matrix_[row][col] = number;
  //   }
  // }
}

void S21Matrix::printMatr() {
  std::cerr << "!!!" << this->matrix_ << "!!!" << std::endl;
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
    S21Matrix pM21 = pM2;  // так работает
    std::cout << pM21.getMatrix(0, 0) << std::endl;
    pM21.printMatr();
    pM2.printMatr();
    S21Matrix matrix2(std::move(pM2));  // вроде робит
    matrix2.printMatr();
    pM2.printMatr();
    pM2.printMatr();
    pM2.printMatr();
    S21Matrix MSum(5, 4);
    MSum.printMatr();
    MSum = matrix2;  // а так не работает (Double free of object 0x7faf23f06350)
    MSum.printMatr();

    // delete pM21;
    // delete pMres;
    // delete pM;
    // delete pM2;
  } catch (const std::exception& err) {
    std::cerr << "[exception] " << err.what() << std::endl;
  }

  return 0;
}