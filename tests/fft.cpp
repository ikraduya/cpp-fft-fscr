#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <algorithm>
#include <limits>
#include <complex>

#include "fft-fscr.hpp"
#include "matrix.hpp"

namespace
{
  using complex_d = std::complex<double>;
  constexpr double absoluteError = 0.00001;

  template<typename T>
  std::vector<T> generate_random_v(size_t size, T lambda = 3.5)
  {
    using value_type = T;
    // We use static in order to instantiate the random engine
    // and the distribution once only.
    // It may provoke some thread-safety issues.

    static std::exponential_distribution<value_type> distribution(lambda);
    static std::default_random_engine generator;

    std::vector<value_type> data(size);
    std::generate(data.begin(), data.end(), []() { return distribution(generator); });
    return data;
  }

  template<typename T, typename U>
  inline void EXPECT_NEAR_V_CPLX(const std::vector<std::complex<T>>& v1, const std::vector<std::complex<U>>& v2) {
    EXPECT_EQ(v1.size(), v2.size());
    for (size_t i = 0; i < v1.size(); ++i) {
      EXPECT_NEAR(v1[i].real(), v2[i].real(), absoluteError);
      EXPECT_NEAR(v1[i].imag(), v2[i].imag(), absoluteError);
    }
  }

  /**
   * Variable used for tests
   */
  static const std::vector<double> input_vector{0.2165, 0.8321, 0.7835, 0.5821, 0.2165, -0.5821, -1.2165, -0.8321};
  static const std::vector<complex_d> expected_result{
    complex_d(0,0), complex_d(0,-3.99998),
    complex_d(0.866,-0.5), complex_d(0,1.91801e-05),
    complex_d(0,0), complex_d(0,-1.91801e-05),
    complex_d(0.866,0.5), complex_d(0,3.99998)
  };

} // anonymous namespace

TEST(DFT, tc1LoopVersion) {
  auto ans = fscr::FFT::dft_loop(input_vector);
  EXPECT_NEAR_V_CPLX(ans, expected_result);
}

TEST(DFT, tc2MatrixVersion) {
  auto ans = fscr::FFT::dft_mat(input_vector);
  EXPECT_NEAR_V_CPLX(ans, expected_result);
}

TEST(DFT, tc3LoopAndMatrixEquality) {
  const std::vector<double> x = generate_random_v<double>(1000);
  EXPECT_NEAR_V_CPLX(fscr::FFT::dft_loop(x), fscr::FFT::dft_mat(x));
}

TEST(FFT, tc1DnC) {
  auto ans = fscr::FFT::fft_dnc(input_vector);
  EXPECT_NEAR_V_CPLX(ans, expected_result);
}

TEST(FFT, tc2DnCInvalid) {
  const std::vector<double> x = generate_random_v<double>(1000);
  EXPECT_THROW(
    fscr::FFT::fft_dnc(x),
    std::invalid_argument);
}

TEST(DFTAndFFT, tc1LoopAndDnCEquality) {
  const std::vector<double> x = generate_random_v<double>(1024);
  EXPECT_NEAR_V_CPLX(fscr::FFT::dft_loop(x), fscr::FFT::fft_dnc(x));
}

TEST(DFTAndFFT, tc2MatrixAndDnCEquality) {
  const std::vector<double> x = generate_random_v<double>(1024);
  EXPECT_NEAR_V_CPLX(fscr::FFT::dft_mat(x), fscr::FFT::fft_dnc(x));
}
