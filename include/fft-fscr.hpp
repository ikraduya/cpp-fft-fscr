#ifndef FSCR_FFT_HPP
#define FSCR_FFT_HPP

#include <vector>
#include <complex>
#include <stdexcept>
#include <type_traits>

#if defined(_WIN32)
#define M_PI   3.14159265358979323846
#else
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include "matrix.hpp"

namespace
{
  /**
   * @brief Slice vector
   */
  template<typename T>
  std::vector<T> segment_vec(const std::vector<T> vec, const size_t start_idx, const size_t finish_idx_one_pass) {
    if (start_idx >= vec.size()) {
      throw std::invalid_argument("segment - start_idx >= vec.size()");
    }
    std::vector<T> res(finish_idx_one_pass - start_idx);
    std::copy(vec.begin() + start_idx, vec.begin() + finish_idx_one_pass, res.begin());
    return res;
  }

  /**
   * @brief Add and multiply vectors values in vec1[i] + (vec2[i] * vec3[i]) manner
   */
  template<typename T>
  std::vector<T> add_multiply_vec(const std::vector<T>& vec1, const std::vector<T>& vec2, const std::vector<T>& vec3) {
    auto res(vec1);
    for (size_t i = 0; i < res.size(); ++i) {
      res[i] += vec2[i] * vec3[i];
    }
    return res;
  }

  /**
   * @brief Concat two vector
   */
  template<typename T>
  std::vector<T> concat_vec(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    auto res(vec1);
    res.insert(res.end(), vec2.begin(), vec2.end());
    return res;
  };

} // anonymous namespace


namespace fscr
{

  class FFT
  {
    public:
    /**
     * @brief Fast Fourier Transform - Divide and Conquer (Cooley-Tukey algorithm)
     */
    template<typename T = double>
    static std::vector<std::complex<T>> fft_dnc(const std::vector<T> &x) {
      static_assert(std::is_floating_point<T>::value, "The input value type should be real");
      
      const size_t N = x.size();
      if (N % 2 != 0) {
        throw std::invalid_argument("Doesn't support vector non power0of-2 size");
      }

      /**
       * @brief Slice odd or even indices from values
       */
      auto slice_step = [](bool isOdd, const std::vector<T> &vec) -> std::vector<T> {
        std::vector<T> res;
        res.reserve(vec.size() / 2);
        
        const size_t start_i = (isOdd) ? 1 : 0;
        for (size_t i = start_i; i < vec.size(); i += 2) {
          res.push_back(vec[i]);
        }
        return res;
      };

      if (N <= 2) {
        return dft_mat(x);
      } else {
        auto X_even = fft_dnc(slice_step(false, x));
        auto X_odd = fft_dnc(slice_step(true, x));
        const auto terms = fscr::expFlatten(
          (std::complex<T>(0, -2*M_PI) * fscr::matrix_d::arange_1d_row(N)) / static_cast<T>(N)
        );
        
        return concat_vec(
          add_multiply_vec(X_even, segment_vec(terms, 0, N / 2), X_odd), 
          add_multiply_vec(X_even, segment_vec(terms, N / 2, N), X_odd) 
        );
      }
    }

    /**
     * @brief Discrete Fourier Transform - Naive Looping
     */
    template<typename T = double>
    static std::vector<std::complex<T>> dft_loop(const std::vector<T> &x) {
      using complex_type = std::complex<T>;
      
      const size_t N = x.size();

      /**
       * @brief DFT for one element
       */
      auto dft_k = [&x, N](const size_t k) -> complex_type {
        complex_type summed(0.0, 0.0);
        for (size_t n=0; n<N; ++n) {
          summed += x[n] * std::exp(complex_type(0, -2) * complex_type(M_PI * k * n) / complex_type(static_cast<T>(N)));
        }
        return summed;
      };

      std::vector<complex_type> res(N);
      for (size_t k=0; k<N; ++k) {
        res[k] = dft_k(k);
      }

      return res;
    }

    /**
     * @brief Discrete Fourier Transform - Matrix operations
     */
    template<typename T = double>
    static std::vector<std::complex<T>> dft_mat(const std::vector<T> &x) {
      const size_t N = x.size();
      const auto n = fscr::matrix_t<T>::arange_1d_row(N);
      const auto k = fscr::matrix_t<T>::arange_1d_col(N);
      constexpr auto min_2_j_pi = std::complex<T>(0, -2 * M_PI);
      const auto M = fscr::exp(((min_2_j_pi * k) * n) / static_cast<T>(N));
      return fscr::dot(M, x);
    }
  };
}



#endif // FSCR_FFT_HPP