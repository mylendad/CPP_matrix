#include <gtest/gtest.h>

#include <cmath>

#include "../s21_matrix_oop.h"

class S21MatrixTest : public ::testing::Test {
 protected:
  void SetUp() override {
    m1_ = S21Matrix(2, 3);
    m1_.SetNumToMatrix(0, 0, 1.0);
    m1_.SetNumToMatrix(0, 1, 2.0);
    m1_.SetNumToMatrix(0, 2, 3.0);
    m1_.SetNumToMatrix(1, 0, 4.0);
    m1_.SetNumToMatrix(1, 1, 5.0);
    m1_.SetNumToMatrix(1, 2, 6.0);

    m2_ = S21Matrix(2, 3);
    m2_.SetNumToMatrix(0, 0, 0.5);
    m2_.SetNumToMatrix(0, 1, 1.5);
    m2_.SetNumToMatrix(0, 2, 2.5);
    m2_.SetNumToMatrix(1, 0, 3.5);
    m2_.SetNumToMatrix(1, 1, 4.5);
    m2_.SetNumToMatrix(1, 2, 5.5);

    square_ = S21Matrix(3, 3);
    square_.SetNumToMatrix(0, 0, 1);
    square_.SetNumToMatrix(0, 1, 2);
    square_.SetNumToMatrix(0, 2, 3);
    square_.SetNumToMatrix(1, 0, 0);
    square_.SetNumToMatrix(1, 1, 4);
    square_.SetNumToMatrix(1, 2, 2);
    square_.SetNumToMatrix(2, 0, 5);
    square_.SetNumToMatrix(2, 1, 2);
    square_.SetNumToMatrix(2, 2, 1);
  }

  S21Matrix m1_;
  S21Matrix m2_;
  S21Matrix square_;
};

TEST_F(S21MatrixTest, DefaultConstructor) {
  S21Matrix m;
  EXPECT_EQ(m.GetRows(), 3);
  EXPECT_EQ(m.GetCols(), 3);
}

TEST_F(S21MatrixTest, ParameterizedConstructor) {
  S21Matrix m(4, 5);
  EXPECT_EQ(m.GetRows(), 4);
  EXPECT_EQ(m.GetCols(), 5);
}

TEST_F(S21MatrixTest, ConstructorNegativeSize) {
  EXPECT_THROW(S21Matrix(-1, 3), std::invalid_argument);
  EXPECT_THROW(S21Matrix(2, -5), std::invalid_argument);
}

TEST_F(S21MatrixTest, CopyConstructor) {
  S21Matrix copy(m1_);
  EXPECT_EQ(copy.GetRows(), 2);
  EXPECT_EQ(copy.GetCols(), 3);
}

TEST_F(S21MatrixTest, MoveConstructor) {
  S21Matrix moved(std::move(m1_));
  EXPECT_EQ(moved.GetRows(), 2);
  EXPECT_EQ(moved.GetCols(), 3);
  EXPECT_EQ(m1_.GetRows(), 0);
  EXPECT_EQ(m1_.GetCols(), 0);
}

TEST_F(S21MatrixTest, AssignmentOperator) {
  S21Matrix copy;
  copy = m1_;
  EXPECT_EQ(copy.GetRows(), 2);
  EXPECT_EQ(copy.GetCols(), 3);
  EXPECT_DOUBLE_EQ(copy.GetMatrix(0, 1), 2.0);
}

TEST_F(S21MatrixTest, MoveAssignmentOperator) {
  S21Matrix moved;
  moved = std::move(m1_);
  EXPECT_EQ(moved.GetRows(), 2);
  EXPECT_EQ(moved.GetCols(), 3);
  EXPECT_EQ(m1_.GetRows(), 0);
}

TEST_F(S21MatrixTest, GetSetElements) {
  EXPECT_DOUBLE_EQ(m1_.GetMatrix(0, 0), 1.0);
  m1_.SetNumToMatrix(1, 1, 99.5);
  EXPECT_DOUBLE_EQ(m1_.GetMatrix(1, 1), 99.5);
}

TEST_F(S21MatrixTest, AccessOperator) {
  m1_(1, 1) = 7.5;
  EXPECT_DOUBLE_EQ(m1_(1, 1), 7.5);
  EXPECT_DOUBLE_EQ(m1_(0, 2), 3.0);
}

TEST_F(S21MatrixTest, AccessOperatorOutOfRange) {
  EXPECT_THROW(m1_(5, 0), std::out_of_range);
  EXPECT_THROW(m1_(0, 5), std::out_of_range);
  EXPECT_THROW(m1_(-1, 0), std::out_of_range);
}

TEST_F(S21MatrixTest, Addition) {
  S21Matrix result = m1_ + m2_;
  EXPECT_DOUBLE_EQ(result(0, 0), 1.5);
  EXPECT_DOUBLE_EQ(result(1, 2), 11.5);
}

TEST_F(S21MatrixTest, AdditionIncompatibleSizes) {
  S21Matrix wrong(3, 3);
  EXPECT_THROW(m1_ + wrong, std::invalid_argument);
}

TEST_F(S21MatrixTest, Subtraction) {
  S21Matrix result = m1_;
  result -= m2_;
  EXPECT_DOUBLE_EQ(result(0, 0), 0.5);
  EXPECT_DOUBLE_EQ(result(1, 2), 0.5);
}

TEST_F(S21MatrixTest, MultiplicationByNumber) {
  m1_.MulNumber(2.5);
  EXPECT_DOUBLE_EQ(m1_(0, 0), 2.5);
  EXPECT_DOUBLE_EQ(m1_(1, 2), 15.0);
}

TEST_F(S21MatrixTest, MatrixMultiplication) {
  S21Matrix a(2, 3);
  a.SetNumToMatrix(0, 0, 1);
  a.SetNumToMatrix(0, 1, 2);
  a.SetNumToMatrix(0, 2, 3);
  a.SetNumToMatrix(1, 0, 4);
  a.SetNumToMatrix(1, 1, 5);
  a.SetNumToMatrix(1, 2, 6);

  S21Matrix b(3, 2);
  b.SetNumToMatrix(0, 0, 7);
  b.SetNumToMatrix(0, 1, 8);
  b.SetNumToMatrix(1, 0, 9);
  b.SetNumToMatrix(1, 1, 10);
  b.SetNumToMatrix(2, 0, 11);
  b.SetNumToMatrix(2, 1, 12);

  // a.MulMatrix(b);
  a *= b;
  EXPECT_EQ(a.GetRows(), 2);
  EXPECT_EQ(a.GetCols(), 2);
  EXPECT_DOUBLE_EQ(a(0, 0), 58);
  EXPECT_DOUBLE_EQ(a(1, 1), 154);
}

TEST_F(S21MatrixTest, MatrixMultiplicationIncompatible) {
  S21Matrix wrong(4, 4);
  EXPECT_THROW(m1_ *= wrong, std::invalid_argument);
}

TEST_F(S21MatrixTest, Equality) {
  S21Matrix copy = m1_;
  EXPECT_TRUE(m1_.EqMatrix(copy) == 0);

  copy.SetNumToMatrix(1, 1, 99.9);
  EXPECT_TRUE(m1_.EqMatrix(copy) == 1);
}

TEST_F(S21MatrixTest, EqualityWithPrecision) {
  S21Matrix a(2, 2);
  a.SetNumToMatrix(0, 0, 0.12345678);
  a.SetNumToMatrix(0, 1, 0.87654321);

  S21Matrix b(2, 2);
  b.SetNumToMatrix(0, 0, 0.12345679);
  b.SetNumToMatrix(0, 1, 0.87654321);

  EXPECT_TRUE(a.EqMatrix(b) == 0);
}

TEST_F(S21MatrixTest, Transpose) {
  S21Matrix transposed = m1_.Transpose();
  EXPECT_EQ(transposed.GetRows(), 3);
  EXPECT_EQ(transposed.GetCols(), 2);
  EXPECT_DOUBLE_EQ(transposed(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(transposed(2, 1), 6.0);
}

TEST_F(S21MatrixTest, Determinant) {
  EXPECT_DOUBLE_EQ(square_.Determinant(), -40.0);

  S21Matrix single(1, 1);
  single(0, 0) = 5.5;
  EXPECT_DOUBLE_EQ(single.Determinant(), 5.5);
}

TEST_F(S21MatrixTest, DeterminantNonSquare) {
  EXPECT_THROW(m1_.Determinant(), std::invalid_argument);
}

TEST_F(S21MatrixTest, CalcComplements) {
  S21Matrix comp = square_.CalcComplements();
  EXPECT_DOUBLE_EQ(comp(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(comp(0, 1), 10.0);
  EXPECT_DOUBLE_EQ(comp(0, 2), -20.0);
}

TEST_F(S21MatrixTest, CalcComplementsNonSquare) {
  EXPECT_THROW(m1_.CalcComplements(), std::invalid_argument);
}

TEST_F(S21MatrixTest, InverseMatrix) {
  S21Matrix inv = square_.InverseMatrix();
  S21Matrix res_square_(3, 3);
  res_square_.SetNumToMatrix(0, 0, 0);
  res_square_.SetNumToMatrix(0, 1, -0.1);
  res_square_.SetNumToMatrix(0, 2, 0.2);
  res_square_.SetNumToMatrix(1, 0, -0.25);
  res_square_.SetNumToMatrix(1, 1, 0.35);
  res_square_.SetNumToMatrix(1, 2, 0.05);
  res_square_.SetNumToMatrix(2, 0, 0.5);
  res_square_.SetNumToMatrix(2, 1, -0.2);
  res_square_.SetNumToMatrix(2, 2, -0.1);

  EXPECT_DOUBLE_EQ(inv(0, 0), res_square_(0, 0));
  EXPECT_DOUBLE_EQ(inv(0, 1), res_square_(0, 1));
  EXPECT_DOUBLE_EQ(inv(0, 2), res_square_(0, 2));
  EXPECT_DOUBLE_EQ(inv(1, 0), res_square_(1, 0));
  EXPECT_DOUBLE_EQ(inv(1, 1), res_square_(1, 1));
  EXPECT_DOUBLE_EQ(inv(1, 2), res_square_(1, 2));
  EXPECT_DOUBLE_EQ(inv(2, 0), res_square_(2, 0));
  EXPECT_DOUBLE_EQ(inv(2, 1), res_square_(2, 1));
  EXPECT_DOUBLE_EQ(inv(2, 2), res_square_(2, 2));
}

TEST_F(S21MatrixTest, InverseSingularMatrix) {
  S21Matrix singular(2, 2);
  singular(0, 0) = 1;
  singular(0, 1) = 2;
  singular(1, 0) = 2;
  singular(1, 1) = 4;
  EXPECT_THROW(singular.InverseMatrix(), std::invalid_argument);
}

TEST_F(S21MatrixTest, CopyToSelf) {
  m1_ = m1_;
  EXPECT_DOUBLE_EQ(m1_(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m1_(1, 2), 6.0);
}

TEST_F(S21MatrixTest, SetRowsCols) {
  m1_.SetRows(3);
  EXPECT_EQ(m1_.GetRows(), 3);
  EXPECT_EQ(m1_.GetCols(), 3);

  m1_.SetCols(4);
  EXPECT_EQ(m1_.GetRows(), 3);
  EXPECT_EQ(m1_.GetCols(), 4);

  EXPECT_DOUBLE_EQ(m1_(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m1_(0, 1), 2.0);
}

TEST_F(S21MatrixTest, LargeMatrixOperations) {
  const int size = 100;
  S21Matrix large(size, size);

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      large(i, j) = i + j;
    }
  }

  S21Matrix transposed = large.Transpose();
  EXPECT_EQ(transposed.GetRows(), size);
  EXPECT_EQ(transposed.GetCols(), size);
  EXPECT_DOUBLE_EQ(transposed(50, 50), large(50, 50));

  large.MulNumber(2.0);
  EXPECT_DOUBLE_EQ(large(99, 99), 396.0);
}

TEST_F(S21MatrixTest, DeterminantPrecision) {
  S21Matrix m(3, 3);
  m(0, 0) = 0.0000001;
  m(0, 1) = 0.0000002;
  m(0, 2) = 0.0000003;
  m(1, 0) = 0.0000004;
  m(1, 1) = 0.0000005;
  m(1, 2) = 0.0000006;
  m(2, 0) = 0.0000007;
  m(2, 1) = 0.0000008;
  m(2, 2) = 0.0000009;

  EXPECT_NEAR(m.Determinant(), 0.0, 1e-15);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
