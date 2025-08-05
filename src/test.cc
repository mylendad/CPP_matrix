#include <math.h>

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

  void Copy(const S21Matrix& other);

  // void Transfer(const S21Matrix&& other);

  void CleanMatrix();

  const S21Matrix& operator=(const S21Matrix& other);

  S21Matrix& operator=(S21Matrix&& other) noexcept;

  void SubMatrix(const S21Matrix& other);

  int EqMatrix(const S21Matrix& other);

  void MulNumber(const double number);

  void MulMatrix(const S21Matrix& other);

  S21Matrix Transpose();

  int GetRows() const;

  int GetCols() const;

  void SetRows(int rows);

  void SetCols(int cols);

  double GetMatrix(int row, int col) const;

  S21Matrix CalcComplements();

  double Determinant();

  S21Matrix InverseMatrix();

  double CalcMinors(S21Matrix& other, int row, int col);

  S21Matrix CreateDeterminateMatrix(S21Matrix& other, int row, int col);

  double EndUnit(S21Matrix& other);

  double TrimerNumb(double src);

  void SetNumToMatrix(int row, int col, double number);

  void PrintMatr();
};

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

S21Matrix::~S21Matrix() { this->CleanMatrix(); }

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  if (this->rows_ != 0 && this->cols_ != 0) {
    this->Copy(other);
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
  // std::cerr << other.matrix_ << "!!!in perenos" << std::endl;
  // this->Transfer(std::move(other));
}

// void S21Matrix::Transfer(const S21Matrix&& other) {
//   rows_ = other.rows_;
//   cols_ = other.cols_;
//   this->matrix_ = other.matrix_;
//   other.matrix_ = nullptr;
//   other.rows_ = 0;
//   other.cols_ = 0;
// }

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::exception();
  } else if (this->matrix_ == nullptr || other.matrix_ == nullptr) {
    throw std::exception();
  } else if (this->rows_ < 1 || this->cols_ < 1 || other.rows_ < 1 ||
             other.cols_ < 1) {
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

const S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (&other == this) return *this;
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->Copy(other);
  return *this;
}

S21Matrix& S21Matrix::operator=(const S21Matrix&& other) noexcept {
  // if (&other == this) return *this;
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = other.matrix_;
  other.matrix_ = nullptr;  // ??
  other.rows_ = 0;          // ??
  other.cols_ = 0;
  return *this;
}

// S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
//     if (this == &other) return *this; // Проверяем самоприсваивание

//     // Освобождаем текущие ресурсы
//     delete[] matrix_; // Освобождаем текущую память

//     // Перемещаем данные из другого объекта
//     rows_ = other.rows_;
//     cols_ = other.cols_;
//     matrix_ = other.matrix_;

//     // Обнуляем перемещаемые ресурсы
//     other.matrix_ = nullptr;
//     other.rows_ = 0;
//     other.cols_ = 0;

//     return *this;
// }

S21Matrix& S21Matrix::operator+(const S21Matrix& other) {
  // creating result matrix

  this->SumMatrix(other);
  // std::cerr << this->matrix_ << "!!!inoperator+" << std::endl;
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->Copy(other);
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  if (this->matrix_ != nullptr) {
    this->CleanMatrix();
  }
  this->Copy(other);
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

void S21Matrix::Copy(const S21Matrix& other) {
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  // std::cerr << this->matrix_ << "!!!Copy begin" << std::endl;
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

  // std::cerr << this->matrix_ << "!!!Copy end" << std::endl;
  // return *this;
}

void S21Matrix::CleanMatrix() {
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
  } else if (this->rows_ < 1 || this->cols_ < 1 || other.rows_ < 1 ||
             other.cols_ < 1) {
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

void S21Matrix::MulNumber(const double number) {
  if (this->matrix_ == nullptr) {
    throw std::exception();
  } else if (this->rows_ < 1 || this->cols_ < 1) {
    throw std::exception();
  } else {
    S21Matrix(this->rows_, this->cols_);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        this->matrix_[i][j] *= number;
      }
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->rows_ != other.cols_ || this->cols_ != other.rows_) {
    throw std::exception();
  } else if (this->matrix_ == nullptr || other.matrix_ == nullptr) {
    throw std::exception();
  } else if (this->rows_ < 1 || this->cols_ < 1 || other.rows_ < 1 ||
             other.cols_ < 1) {
    throw std::exception();
  } else {
    S21Matrix(this->rows_, this->cols_);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        for (int k = 0; k < other.rows_; k++)
          this->matrix_[i][j] += (this->matrix_[i][k] * other.matrix_[k][j]);
      }
    }
  }
}

S21Matrix S21Matrix::Transpose() {
  if (this->matrix_ == nullptr) {
    throw std::exception();
  } else if (this->rows_ < 1 || this->cols_ < 1) {
    throw std::exception();
  } else {
    S21Matrix res(this->cols_, this->rows_);
    for (int i = 0; i < this->cols_; i++) {
      for (int j = 0; j < this->rows_; j++) {
        res.matrix_[i][j] = this->matrix_[j][i];
      }
    }
    return res;
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
      if (TrimerNumb(this->matrix_[i][j]) == TrimerNumb(other.matrix_[i][j])) {
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

int S21Matrix::GetRows() const { return this->rows_; }

int S21Matrix::GetCols() const { return this->cols_; }

void S21Matrix::SetRows(int rows) {
  this->rows_ = rows;
  S21Matrix(this->rows_, this->cols_);
}

void S21Matrix::SetCols(int cols) {
  this->rows_ = cols;
  S21Matrix(this->rows_, this->cols_);
}

double S21Matrix::GetMatrix(int row, int col) const {
  return this->matrix_[row][col];
}

void S21Matrix::SetNumToMatrix(int row, int col, double number) {
  this->matrix_[row][col] = number;
}

void S21Matrix::PrintMatr() {
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

S21Matrix S21Matrix::CalcComplements() {
  if (this->matrix_ == nullptr) {
    throw std::exception();
  } else if (this->rows_ != this->cols_) {
    throw std::exception();
  } else if (this->rows_ <= 1) {
    throw std::exception();
  }
  // int res = 0;
  // if (this->rows_ == 1) res = 1;
  // int flag = 0;
  S21Matrix res(this->rows_);
  std::cerr << "!!! пока все ок!;" << std::endl;
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
        std::cerr << "!!! j = " << j << std::endl;
        res.matrix_[i][j] =

            CalcMinors(*this, i, j) *
            pow(-1, ((this->rows_ + 1 - i) + (this->rows_ + 1 - j)));
      }
    }
    // res = 0;
  }
  return res;
}

S21Matrix S21Matrix::CreateDeterminateMatrix(S21Matrix& other, int row,
                                             int col) {
  // std::cerr << "!!! begin CreateDeterminateMatrix(;" << std::endl;
  // std::cerr << "!!!rows = ;" << rows << std::endl;
  if (other.matrix_ == nullptr) {
    throw std::exception();

  } else if (other.rows_ < 1) {
    throw std::exception();
  }
  S21Matrix res(other.rows_ - 1);
  std::cout << "!!! other.rows_" << other.rows_ << std::endl;
  // std::cerr << "!!! after create res;" << std::endl;
  for (int i = 0, x = 0; i < other.rows_; i++) {  // change x, y
    for (int j = 0, y = 0; j < other.rows_; j++) {
      if (j != col && i != row) {
        res.matrix_[x][y] = other.matrix_[i][j];
        // std::cerr << "!!! y == " << y << std::endl;  // change x, y
        y++;
      }
      if (y == other.rows_ - 1) {
        y = 0;
        x++;
      }
    }
  }
  // std::cerr << "!!! end" << std::endl;
  return res;
}

double S21Matrix::CalcMinors(S21Matrix& other, int row, int col) {
  std::cout << "!!!other.matrix_ = " << other.matrix_ << std::endl;

  if (other.matrix_ == nullptr) {
    throw std::exception();

  } else if (this->rows_ < 1) {
    throw std::exception();
  }
  // std::cout << "!!! начало CalcMinors!;" << std::endl;
  // std::cout << "!!! int rows= " << row << std::endl;
  double minor = 0.0;
  int size_matrix = this->rows_;
  // std::cerr << "!!! size_matrix = this->rows_= " << this->rows_ << std::endl;
  // S21Matrix Copy{other};
  S21Matrix temp;
  // std::cout << "!!! все еще ок;" << std::endl;
  if (other.rows_ > 2) {
    temp = CreateDeterminateMatrix(other, row, col);
    for (int j = 0; j < size_matrix - 1; j++) {
      double a = CalcMinors(temp, 0, j) * temp.matrix_[0][j];
      if (j % 2 > 0) a *= -1;
      minor += a;
    }
    // std::cout << "!!! и щас?" << std::endl;
  } else if (other.rows_ == 2)
    std::cout << "!!! other.rows_ == 2" << other.rows_ << std::endl;
  minor = EndUnit(temp);

  // delete temp;
  std::cout << "!!! end calc minors" << other.rows_ << std::endl;
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
  char buffer[20];  // add define

  sprintf(buffer, "%.7lf", src);
  sscanf(buffer, "%lf", &result);
  return result;
}

double S21Matrix::Determinant() {
  if (this->rows_ != this->cols_) {
    throw std::exception();
  } else if (this->matrix_ == nullptr) {
    throw std::exception();
  } else if (this->rows_ < 1) {
    throw std::exception();
  }

  double result = 0;
  if (this->rows_ == 1)
    result = this->matrix_[0][0];
  else {
    if (this->rows_ == 2)
      result = EndUnit(*this);
    else {
      for (int y = 0; y < this->rows_; y++) {
        double a = CalcMinors(*this, 0, y) * this->matrix_[0][y];
        if (y % 2 > 0) a *= -1;
        result += a;
      }
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (this->rows_ != this->cols_) {
    throw std::exception();
  } else if (this->matrix_ == nullptr) {
    throw std::exception();
  } else if (this->rows_ < 1) {
    throw std::exception();
  }
  double determin = this->Determinant();
  if (determin == 0) throw std::exception();
  // s21_create_matrix(A->rows, A->columns, result);
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

int main() {
  try

  {
    std::cerr << "!!!start" << std::endl;
    S21Matrix pM(2, 2);

    S21Matrix pM2(2, 2);  //
    // S21Matrix* pMres = new S21Matrix(4, 5);
    // pM->SumMatrix->matrix_[0][0] = 1.0;
    pM2.SetNumToMatrix(0, 0, 1.0);
    pM2.PrintMatr();
    pM.PrintMatr();
    pM.SumMatrix(pM2);
    std::cout << pM.GetMatrix(0, 0) << std::endl;
    S21Matrix pMres{pM};  // Copy
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
    pM22 = MSum.Transpose();
    std::cout << "!!! MSum.Transpose();" << std::endl;
    std::cout << "!!!pM22" << std::endl;
    pM22.PrintMatr();
    std::cout << "!!! MSum.GetRows() " << MSum.GetRows() << std::endl;
    // S21Matrix pM23 = MSum.CalcComplements();
    // std::cerr << "!!!pM23 = MSum.CalcComplements();" << std::endl;
    // std::cerr << "!!!pM23" << std::endl;
    // pM23.PrintMatr();

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