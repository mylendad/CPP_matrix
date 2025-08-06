#include <math.h>

#include <cstring>
#include <iostream>
#include <stdexcept>

#define TRUE 1
#define FALSE 0
#define BUFFER_TRIMMER 20

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

  void SubMatrix(const S21Matrix& other);

  int EqMatrix(const S21Matrix& other);

  void MulNumber(const double number);

  void MulMatrix(const S21Matrix& other);

  S21Matrix Transpose();

  S21Matrix CalcComplements();

  double Determinant();

  S21Matrix InverseMatrix();

  const S21Matrix& operator=(const S21Matrix& other);

  S21Matrix& operator=(S21Matrix&& other) noexcept;

  S21Matrix& operator+(const S21Matrix& other);

  S21Matrix& operator+=(const S21Matrix& other);

  S21Matrix& operator-=(const S21Matrix& other);

  S21Matrix& operator*=(const S21Matrix& other);

  double& operator()(int i, int j);

  int GetRows() const;

  int GetCols() const;

  void SetRows(int rows);

  void SetCols(int cols);

  double GetMatrix(int row, int col) const;

  void CopyMatrix(const S21Matrix& other);

  void CleanMatrix();

  double CalcMinors(S21Matrix& other, int row, int col);

  S21Matrix CreateDeterminateMatrix(S21Matrix& other, int row, int col);

  double EndUnit(S21Matrix& other);

  double TrimerNumb(double src);

  bool IsEvenNumber(int number);

  void SetNumToMatrix(int row, int col, double number);

  void PrintMatr();
};