#include <cstring>
#include <iostream>
#include <stdexcept>

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_ = nullptr;

 public:
  // S21Matrix() : rows_{3}, cols_{3} {
  //   this->matrix_ = new double*[rows_];
  //   for (int i = 0; i < this->rows_; i++) {
  //     matrix_[i] = new double[this->cols_]{0};
  //   }
  // }

  // S21Matrix(int rows) : rows_{rows}, cols_{rows} {
  //   // this->rows_ = rows;
  //   // this->cols_ = rows;
  //   // this->matrix_ = new double*[rows * rows]();
  //   this->matrix_ = new double*[rows];
  //   for (int i = 0; i < rows; i++) {
  //     matrix_[i] = new double[rows]{0};
  //   }
  // }

  S21Matrix(int rows, int cols) : rows_{rows}, cols_{cols} {
    this->matrix_ = new double*[rows];
    for (int i = 0; i < rows; i++) {
      this->matrix_[i] = new double[cols]{0};
    }
  }

  S21Matrix(int rows) : S21Matrix(rows, rows) {}

  S21Matrix() : S21Matrix(3, 3) {}

  ~S21Matrix() { delete[] matrix_; }

  S21Matrix(const S21Matrix& other) : rows_(other.rows_), cols_(other.cols_) {
    // if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    //   throw std::exception();
    // }
    std::cerr << "hi" << std::endl;
    this->matrix_ = new double*[other.rows_];  // вынести
    for (int i = 0; i < other.rows_; i++) {
      matrix_[i] = new double[this->cols_]{0};
      std::memcpy(this->matrix_[i], other.matrix_[i],
                  other.cols_ * sizeof(double));
      // std::cerr << "!" << *this->matrix_[i] << "!" << " !" <<
      // *other.matrix_[i]
      // << "!" << std::endl;
    }
    std::cerr << this->matrix_ << "!!!" << std::endl;
    std::cerr << other.matrix_ << "!!!" << std::endl;
  }

  S21Matrix(const S21Matrix&& other) : rows_(other.rows_), cols_(other.cols_) {
    // if (this->rows_ != other.rows_ && this->cols_ != other.cols_)
    this->matrix_ = new double*[other.rows_];  // вынести
    for (int i = 0; i < other.rows_; i++) {
      matrix_[i] = new double[this->cols_]{0};
      std::memcpy(this->matrix_[i], other.matrix_[i],
                  other.cols_ * sizeof(double));
      // std::cerr << "!" << *this->matrix_[i] << "!" << " !" <<
      // *other.matrix_[i]
      // << "!" << std::endl;
    }

    // } else {
    this->matrix_ = new double*[other.rows_];
    for (int i = 0; i < other.rows_; i++) {
      matrix_[i] = new double[this->cols_]{0};
      std::memcpy(this->matrix_[i], other.matrix_[i],
                  other.cols_ * sizeof(double));
    }
    // }
    for (int i = 0; i < this->rows_; i++) {
      delete[] other.matrix_[i];
    }
  }

  void SumMatrix(const S21Matrix& other) {
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

  void operator+(const S21Matrix& other) {
    // creating result matrix

    this->SumMatrix(other);
  }

  const S21Matrix& operator=(
      const S21Matrix& left) {  // спросить как это работает
    if (&left == this) return *this;
    if (this->matrix_ != nullptr) {
      for (int i = 0; i < this->rows_; i++) {
        delete[] this->matrix_[i];
      }
      delete[] this->matrix_;
    }
    *this = left;  // испраить или оставить так?
    return *this;
  }

  void SubMatrix(const S21Matrix& other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_ ||
        this->matrix_ == nullptr || other.matrix_ == nullptr) {
      throw std::exception();

    }
    // else if (this->matrix_ == nullptr || other.matrix_ == nullptr) {
    //   throw std::exception();
    // }
    else {
      S21Matrix(this->rows_, this->cols_);
      for (int i = 0; i < this->rows_; i++) {
        for (int j = 0; j < this->cols_; j++) {
          this->matrix_[i][j] -= other.matrix_[i][j];
        }
      }
    }
  }

  void EqMatrix(const S21Matrix& other) {
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

  int getRows() const { return this->rows_; }

  int getCols() const { return this->cols_; }

  void setRows(int rows) { this->rows_ = rows; }

  void setCols(int cols) { this->cols_ = cols; }

  // double** setNumToMatrix(int rows, int cols, double number,
  //                    const S21Matrix& other) {
  //   for (int i = 0; i < rows; i++) {
  //     for (int j = 0; j < cols; j++) {
  //       other.matrix_[i][j] = number;
  //     }
  //   }
  // }

  double getMatrix(int row, int col) const {
    // for (int i = 0; i < this->rows_; i++) {
    //   for (int j = 0; j < this->cols_; j++) {
    return this->matrix_[row][col];
    //   }
    // }
  }

  void setNumToMatrix(int row, int col, double number) {
    // for (int i = 0; i < this->rows_; i++) {
    //   for (int j = 0; j < this->cols_; j++) {
    this->matrix_[row][col] = number;
    //   }
    // }
  }

  void printMatr() {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        if (i == this->rows_ - 1 && j == this->cols_ - 1)
          std::cout << this->matrix_[i][j];
        else if (j == this->cols_ - 1)
          std::cout << this->matrix_[i][j] << std::endl;
        else
          std::cout << this->matrix_[i][j] << " ";
      }
    }
  }
};

int main() {
  try

  {
    S21Matrix* pM = new S21Matrix(5, 4);

    S21Matrix* pM2 = new S21Matrix(5, 4);  //
    // S21Matrix* pMres = new S21Matrix(4, 5);
    // pM->SumMatrix->matrix_[0][0] = 1.0;
    pM2->setNumToMatrix(0, 0, 1.0);
    pM->SumMatrix(*pM2);
    std::cout << pM->getMatrix(0, 0) << std::endl;
    S21Matrix pMres{*pM};  // copy
    pMres.printMatr();
    std::cout << pMres.getMatrix(0, 0) << std::endl;
    pM->SubMatrix(*pM2);
    std::cout << pM->getMatrix(0, 0) << std::endl;
    S21Matrix* pM21 = pM2;
    std::cout << pM21->getMatrix(0, 0) << std::endl;
    ы pM21->printMatr();
    S21Matrix matrix2(std::move(*pM2));  // вроде робит
    matrix2.printMatr();

    // delete pM21;
    // delete pMres;
    delete pM;
    delete pM2;
  } catch (const std::exception& err) {
    std::cerr << "[exception] " << err.what() << std::endl;
  }

  return 0;
}