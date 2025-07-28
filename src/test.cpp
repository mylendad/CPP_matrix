#include <iostream>
#include <stdexcept>

class S21Matrix {
 private:
  int rows_{0}, cols_{0};
  double** matrix_ = NULL;

 public:
  S21Matrix() {
    this->rows_ = 3;
    this->cols_ = 3;
    this->matrix_ = new double*[this->rows_ * this->cols_]();
    for (int i = 0; i < this->rows_; i++) {
      matrix_[i] = new double[this->cols_]{0};
    }
  };

  S21Matrix(int rows) {
    this->rows_ = rows;
    this->cols_ = rows;
    // this->matrix_ = new double*[rows * rows]();
    this->matrix_ = new double*[rows];
    for (int i = 0; i < rows; i++) {
      matrix_[i] = new double[rows]{0};
    }
  };

  S21Matrix(int rows, int cols) {
    this->rows_ = rows;
    this->cols_ = cols;
    // this->matrix_ = new double*[rows * rows]();
    this->matrix_ = new double*[rows];
    for (int i = 0; i < rows; i++) {
      matrix_[i] = new double[cols]{0};
    }
  };

  ~S21Matrix() { delete[] matrix_; }

  void SumMatrix(const S21Matrix& other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
      throw std::exception();

    } else if (this->matrix_ == NULL || other.matrix_ == NULL) {
      throw std::exception();
    }
    S21Matrix(this->rows_, this->cols_);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        this->matrix_[i][j] += other.matrix_[i][j];
      }
    }
  }

  // S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  //   // creating result matrix
  //   S21Matrix res(this->rows_, this->cols_);
  //   res.SumMatrix(other);
  //   return res;
  // }

  void SubMatrix(const S21Matrix& other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
      throw std::exception();

    } else if (this->matrix_ == NULL || other.matrix_ == NULL) {
      throw std::exception();
    }
    S21Matrix(this->rows_, this->cols_);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        this->matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  }

  void EqMatrix(const S21Matrix& other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
      throw std::exception();

    } else if (this->matrix_ == NULL || other.matrix_ == NULL) {
      throw std::exception();
    }
    S21Matrix(this->rows_, this->cols_);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        this->matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  }

  int getRows() { return this->rows_; }

  int getCols() { return this->cols_; }

  void setRows(int rows) { this->rows_ = rows; }

  void setCols(int cols) { this->cols_ = cols; }

  // double** setMatrix(int rows, int cols, double number,
  //                    const S21Matrix& other) {
  //   for (int i = 0; i < rows; i++) {
  //     for (int j = 0; j < cols; j++) {
  //       other.matrix_[i][j] = number;
  //     }
  //   }
  // }

  double getMatrix(int row, int col) {
    // for (int i = 0; i < this->rows_; i++) {
    //   for (int j = 0; j < this->cols_; j++) {
    return this->matrix_[row][col];
    //   }
    // }
  }

  void setMatrix(int row, int col, double number) {
    // for (int i = 0; i < this->rows_; i++) {
    //   for (int j = 0; j < this->cols_; j++) {
    this->matrix_[row][col] = number;
    //   }
    // }
  }

  void printMatr() {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        std::cout << this->matrix_[i][j] << std::endl;
      }
    }
  }
};

int main() {
  try

  {
    S21Matrix* pM = new S21Matrix(4, 5);
    S21Matrix* pM2 = new S21Matrix(4);  //
    // S21Matrix* pMres = new S21Matrix(4, 5);
    // pM->SumMatrix->matrix_[0][0] = 1.0;
    pM2->setMatrix(0, 0, 1);
    pM->SumMatrix(*pM2);
    std::cout << pM->getMatrix(0, 0) << std::endl;
    pM->SubMatrix(*pM2);
    std::cout << pM->getMatrix(0, 0) << std::endl;
    // pM->printMatr();
    delete pM;
    delete pM2;
  } catch (const std::exception& err) {
    std::cerr << "[exception] " << err.what() << std::endl;
  }
  return 0;

  // delete pMres;
}