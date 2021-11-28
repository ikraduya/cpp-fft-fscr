#include <gtest/gtest.h>
#include <stdexcept>

#include "matrix.hpp"

namespace
{
  template<typename T>
  std::vector<T> generate_random(size_t size)
  {
    using value_type = T;
    // We use static in order to instantiate the random engine
    // and the distribution once only.
    // It may provoke some thread-safety issues.

    static std::uniform_real_distribution<value_type> distribution(
        std::numeric_limits<value_type>::min(),
        std::numeric_limits<value_type>::max());
    static std::default_random_engine generator;

    std::vector<value_type> data(size);
    std::generate(data.begin(), data.end(), []() { return distribution(generator); });
    return data;
  }

  using complex_d = std::complex<double>;
  using complex_f = std::complex<float>;
}

TEST(MatrixConstructor, tc1DefaultType) {
  fscr::matrix mat(2, 3);
  EXPECT_EQ(mat.rows, 2);
  EXPECT_EQ(mat.cols, 3);
}

TEST(MatrixConstructor, tc2DoubleType) {
  fscr::matrix_d mat({{1, 2, 3}, {1, 2, 3}});
  EXPECT_EQ(mat.rows, 2);
  EXPECT_EQ(mat.cols, 3);
}

TEST(MatrixConstructor, tc3FloatType) {
  fscr::matrix_f mat({{1, 2, 3}, {1, 2, 3}});
  EXPECT_EQ(mat.rows, 2);
  EXPECT_EQ(mat.cols, 3);
}
TEST(MatrixConstructor, tc4ComplexDoubleType) {
  fscr::matrix_cd mat({
    {complex_d(1, 1), complex_d(2, 2), complex_d(3, 3)},
    {complex_d(1, -1), complex_d(2, -2), complex_d(3, -3)}
  });
  EXPECT_EQ(mat.rows, 2);
  EXPECT_EQ(mat.cols, 3);
}
TEST(MatrixConstructor, tc5ComplexFloatType) {
  fscr::matrix_cd mat({
    {complex_f(1, 1), complex_f(2, 2), complex_f(3, 3)},
    {complex_f(1, -1), complex_f(2, -2), complex_f(3, -3)}
  });
  EXPECT_EQ(mat.rows, 2);
  EXPECT_EQ(mat.cols, 3);
}

TEST(MatrixConstructor, tc6InvalidRows) {
  EXPECT_THROW(
    fscr::matrix_cd mat(0, 1),
    std::invalid_argument);
}

TEST(MatrixConstructor, tc7InvalidCols) {
  EXPECT_THROW(
    fscr::matrix_cd mat(1, 0),
    std::invalid_argument);
}

TEST(MatrixConstructor, tc8NonSquareMatrix) {
  EXPECT_THROW(
    fscr::matrix_cd mat({
      {1, 2, 3}, {3, 4}
    }),
    std::invalid_argument);
}

TEST(MatrixRandom, tc1InvalidRows) {
  EXPECT_THROW(
    fscr::matrix_d::random(0, 3),
    std::invalid_argument);
}

TEST(MatrixRandom, tc2InvalidCols) {
  EXPECT_THROW(
    fscr::matrix_d::random(2, 0),
    std::invalid_argument);
}

TEST(MatrixRandom, tc3DoubleType) {
  fscr::matrix_d mat = fscr::matrix_d::random(2, 3);
  EXPECT_EQ(mat.rows, 2);
  EXPECT_EQ(mat.cols, 3);
}

TEST(MatrixRandom, tc4FloatType) {
  fscr::matrix_f mat = fscr::matrix_f::random(2, 3);
  EXPECT_EQ(mat.rows, 2);
  EXPECT_EQ(mat.cols, 3);
}

TEST(MatrixSubscript, tc1SquareBracket) {
  fscr::matrix_d mat({{1, 2, 3}, {4, 5, 6}});
  EXPECT_EQ(mat[0][0], 1);
  EXPECT_EQ(mat[0][1], 2);
  EXPECT_EQ(mat[0][2], 3);
  EXPECT_EQ(mat[1][0], 4);
  EXPECT_EQ(mat[1][1], 5);
  EXPECT_EQ(mat[1][2], 6);
}

TEST(MatrixSubscript, tc2RoundBracket) {
  fscr::matrix_d mat({{1, 2, 3}, {4, 5, 6}});
  EXPECT_EQ(mat(0, 0), 1);
  EXPECT_EQ(mat(0, 1), 2);
  EXPECT_EQ(mat(0, 2), 3);
  EXPECT_EQ(mat(1, 0), 4);
  EXPECT_EQ(mat(1, 1), 5);
  EXPECT_EQ(mat(1, 2), 6);
}

TEST(MatrixMultiplication, tc1SameSize) {
  fscr::matrix_d mat1({{1, 1, 1}, {2, 2, 2}, {3, 3, 3}});
  fscr::matrix_d mat2({{1, 1, 1}, {2, 2, 2}, {3, 3, 3}});
  auto result = mat1 * mat2;
  fscr::matrix_d expected_result({{6, 6, 6}, {12, 12, 12}, {18, 18, 18}});
  EXPECT_EQ(result, expected_result);
}

TEST(MatrixMultiplication, tc2NonSameSizeValid) {
  fscr::matrix_d mat1({{1, 2}, {2, 4}, {3, 5}});
  fscr::matrix_d mat2({{1,}, {2,}});
  auto result = mat1 * mat2;
  fscr::matrix_d expected_result({{5,}, {10,}, {13,}});
  EXPECT_EQ(result, expected_result);
}

TEST(MatrixMultiplication, tc3InvalidSizes) {
  fscr::matrix_d mat1({{1, 2, 3}, {2, 4, 5}, {3, 5, 6}});
  fscr::matrix_d mat2({{1,}, {2,}});
  EXPECT_THROW(
    mat1 * mat2,
    std::invalid_argument);
}

TEST(MatrixMultiplication, tc4ComplexTimesMatrix) {
  complex_d cplx(10, -1);
  fscr::matrix_d mat1({{1, 2}, {2, 4}, {3, 5}});
  auto result = cplx * mat1;
  fscr::matrix_cd expected_result({
    {complex_d(10, -1), complex_d(20, -2)}, 
    {complex_d(20, -2), complex_d(40, -4)}, 
    {complex_d(30, -3), complex_d(50, -5)}
  });
  EXPECT_EQ(result, expected_result);
}
