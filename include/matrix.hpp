#ifndef FSCR_MATRIX_HPP
#define FSCR_MATRIX_HPP

#include <initializer_list>
#include <type_traits>
#include <complex>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>

namespace
{
  /**
   * Check complex type
   */
  template<typename T>
  struct is_complex_t : public std::false_type {};
  
  template<typename T>
  struct is_complex_t<std::complex<T>> : public std::true_type {};
  
  template<typename T>
  constexpr bool is_complex() { return is_complex_t<T>::value; }

  /**
   * @brief Check sub-initializer lists size should be equal
   */
  template<typename T>
  bool isSubInitializerListsEqualSize(const std::initializer_list<std::initializer_list<T>> &lis) {
    auto sub_list_itr = lis.begin();
    const size_t sub_list_size = sub_list_itr->size();

    bool is_equal_size = true;
    ++sub_list_itr;
    while (sub_list_itr != lis.end())
    {
      if (sub_list_itr->size() != sub_list_size) {
        is_equal_size = false;
        break;
      }
      ++sub_list_itr;
    }
    return is_equal_size;
  }
}

namespace fscr
{ 
  using complex_d = std::complex<double>;
  
  template<typename T = double>
  class matrix_t {
    private:
    std::vector<T> mem;
    
    public:
    size_t rows = 0, cols = 0;

    /**
     * Section - Constructors
     */
    matrix_t() = default;
    
    matrix_t(const size_t rows, const size_t cols) : rows{rows}, cols{cols} {
      if (rows == 0) {
        throw std::invalid_argument("fscr::matrix_t::matrix_t - rows must be > 0");
      } else if (cols == 0) {
        throw std::invalid_argument("fscr::matrix_t::matrix_t - cols must be > 0");
      }
      mem = std::vector<T>(rows * cols);
    }
    
    matrix_t(const size_t rows, const size_t cols, const T default_val) : rows{rows}, cols{cols} {
      if (rows == 0) {
        throw std::invalid_argument("fscr::matrix_t::matrix_t - rows must be > 0");
      } else if (cols == 0) {
        throw std::invalid_argument("fscr::matrix_t::matrix_t - cols must be > 0");
      }
      mem = std::vector<T>(rows * cols, default_val);
    }
    
    matrix_t(const std::initializer_list<std::initializer_list<T>> &mat) {
      static_assert(std::is_arithmetic<T>::value || is_complex<T>(), "The value should be number (integer, real, or complex)");
      
      if (mat.size() > 0) {
        if (!isSubInitializerListsEqualSize(mat)) {
          throw std::invalid_argument("fscr::matrix_t - Each rows should be the same size");
        }
        
        this->rows = mat.size();
        this->cols = mat.begin()->size();

        mem.reserve(rows * cols);
        for (const auto& m: mat) {
          mem.insert(mem.end(), m.begin(), m.end());
        }
      }
    }

    /**
     * Section - Move and copy constructor
     */
    #if defined(_WIN32) && (_MSC_VER <= 1900)
    // Visual Studio 2015 and earlier do not support "= default" argument
    matrix_t(matrix_t&& source) : 
      mem{ std::move(source.mem) }, rows{ source.rows }, cols{ source.cols } {}
    
    matrix_t(const matrix_t& source) : 
      mem{ source.mem }, rows{ source.rows }, cols{ source.cols } {}
    #else
    matrix_t(matrix_t&& source) = default;
    
    matrix_t(const matrix_t& source) = default;
    #endif

    /**
     * Section - Operator overloadings
     */
    bool operator==(const matrix_t<T>& other) const {
      if (this->rows != other.rows || this->cols != other.cols) {
        return false;
      }
      return (this->mem == other.mem);
    }

    /**
     * @brief Selector to allow mat[i][j] style
     */
    T* operator[](const size_t i_row) {
      return mem.data() + (i_row * this->cols);
    }

    /**
     * @brief Selector to allow mat(i, j) style
     */
    T& operator()(const size_t i_row, const size_t i_col) {
      return mem[i_row * this->cols + i_col];
    }

    /**
     * @brief Matrix multiplication
     */
    matrix_t<T> operator*(const matrix_t<T>& other) const {
      if (this->cols != other.rows) {
        throw std::invalid_argument("fscr::matrix_t::operator* - Matrices dimensions does not meet the valid condition for matrix multiplication");
      }
      
      matrix_t<T> res(this->rows,  other.cols, 0);
      for (size_t i = 0; i < this->rows; ++i) {
        for (size_t j = 0; j < other.cols; ++j) {
          for (size_t k = 0; k < this->cols; ++k) {
            res[i][j] += this->mem[i * this->cols + k] * other.mem[k * other.cols + j];
          }
        }
      }

      return res;
    }

    /**
     * Section - Slice function
     */
    matrix_t<T> block(const size_t start_row, const size_t start_col, const size_t row_size, const size_t col_size) const {
      if (start_row >= this->rows) {
        throw std::invalid_argument("fscr::matrix_t::block - start_row >= matrix rows");
      } else if (start_col >= this->col) {
        throw std::invalid_argument("fscr::matrix_t::block - start_col >= matrix cols");
      }
      if (start_row + row_size >= this->rows) {
        throw std::invalid_argument("fscr::matrix_t::block - start_row + row_size >= matrix rows");
      } else if (start_col + col_size >= this->cols) {
        throw std::invalid_argument("fscr::matrix_t::block - start_row + col_size >= matrix cols");
      }

      matrix_t<T> res(row_size, col_size);
      for (size_t i = 0; i < row_size; ++i) {
        auto start_idx = (start_row + i) * this->cols + start_col;
        std::copy(
          this->mem.begin() + start_idx,
          this->mem.begin() + (start_idx + col_size),
          res.mem.begin() + (i * col_size)
        );
      }

      return res;
    }

    /**
     * Section - Friends Operator overloadings
     */
    template<typename S, typename U>
    friend matrix_t<complex_d> operator*(const std::complex<U>& elem, const matrix_t<S>& mat);
    
    template<typename S, typename U>
    friend matrix_t<complex_d> operator*(const matrix_t<std::complex<S>>& elem, const matrix_t<U>& mat);

    template<typename S>
    friend matrix_t<complex_d> operator/(const matrix_t<std::complex<S>>& mat, const double elem);

    /**
     * Section - Generator function
     */

    /**
     * @brief Produce random value matrix with uniform distribution
     */
    static matrix_t<T> random(const size_t rows, const size_t cols) {
      static_assert(!is_complex<T>(), "Currently fscr::matrix_t::random can't create random complex-value matrix");
      if (rows == 0 || cols == 0) {
        throw std::invalid_argument("fscr::matrix_t::random - row and col size should be > 0");
      }
      matrix_t<T> res(rows, cols);
      
      static std::uniform_real_distribution<T> distribution(
        std::numeric_limits<T>::min(),
        std::numeric_limits<T>::max());
      static std::default_random_engine generator;

      std::generate(res.mem.begin(), res.mem.end(), []() {
        return distribution(generator);
      });

      return std::move(res);
    }

    /**
     * @brief Produce 1 row matrix with evenly spaced value [0, stop)
     * 
     * @param stop end of value (does not include this value)
     */
    static matrix_t<T> arange_1d_row(const size_t stop) {
      matrix_t<T> res(1, stop);

      auto mem_itr = res.mem.begin();
      for (size_t i = 0; i < stop; ++i) {
        *mem_itr = static_cast<T>(i);
        ++mem_itr;
      }

      return res;
    }

    /**
     * @brief Produce 1 col matrix with evenly spaced value [0, stop)
     * 
     * @param stop end of value (does not include this value)
     */
    static matrix_t<T> arange_1d_col(const size_t stop) {
      matrix_t<T> res(stop, 1);

      auto mem_itr = res.mem.begin();
      for (size_t i = 0; i < stop; ++i) {
        *mem_itr = static_cast<T>(i);
        ++mem_itr;
      }

      return res;
    }
    

    /**
     * Section - Friends Math operations
     */
    /**
     * @brief Dot operation, dot(matrix, vector)
     * @return std::vector<complex_d> return vector of complex<double>
     */
    template<typename S, typename U>
    friend std::vector<complex_d> dot(const matrix_t<S>& mat, const std::vector<U>& vec);
    
    /**
     * @brief Apply base-2 exponent to each matrix value
     */
    template<typename S>
    friend matrix_t<S> exp(const matrix_t<S>& mat);

    /**
     * @brief Apply base-2 exponent to each matrix value and convert it to vector
     * 
     * @param mat must be 1-D matrix (row or col)
     */
    template<typename S>
    friend std::vector<complex_d> expFlatten(const matrix_t<std::complex<S>> &mat);

    /**
     * @brief Printout matrix values
     */
    void print() const {
      std::cout << '[';
      for (size_t r = 0; r < rows; ++r) {
        if (r > 0) {
          std::cout << ' ';
        }
        for (size_t c = 0; c < cols; ++c) {
          std::cout << mem[r * cols + c];
          if (r != rows-1 || c != cols-1) {
            std::cout << ',';
            if (c < cols-1) {
              std::cout << ' ';
            }
          }
        }
        if (r < rows-1) {
          std::cout << '\n';
        }
      }
      std::cout << ']' << std::endl;
    }
  };

  /**
   * Section - Friends functions implementation
   */
  template<typename S, typename U>
  matrix_t<complex_d> operator*(const std::complex<U>& elem, const matrix_t<S>& mat) {
    matrix_t<complex_d> res(mat.rows, mat.cols);
    for (size_t i = 0; i < mat.mem.size(); ++i) {
      res.mem[i] = mat.mem[i] * elem;
    }
    
    return std::move(res);
  }

  template<typename S, typename U>
  matrix_t<complex_d> operator*(const matrix_t<std::complex<S>>& matc, const matrix_t<U>& mat) {
    if (matc.cols != mat.rows) {
      throw std::invalid_argument("fscr::operator* - Matrices dimensions does not meet the valid condition for matrix multiplication");
    }
    
    matrix_t<complex_d> res(matc.rows,  mat.cols, complex_d());
    for (size_t i = 0; i < matc.rows; ++i) {
      for (size_t j = 0; j < mat.cols; ++j) {
        for (size_t k = 0; k < matc.cols; ++k) {
          res[i][j] += matc.mem[i * matc.cols + k] * mat.mem[k * mat.cols + j];
        }
      }
    }

    return res;
  }

  template<typename S>
  matrix_t<complex_d> operator/(const matrix_t<std::complex<S>>& mat, const double elem) {
    matrix_t<complex_d> res(mat.rows, mat.cols);
    for (size_t i = 0; i < mat.mem.size(); ++i) {
      res.mem[i] = mat.mem[i] / elem;
    }

    return res;
  }

  template<typename S, typename U>
  std::vector<complex_d> dot(const matrix_t<S>& mat, const std::vector<U>& vec) {
    if (mat.cols != vec.size()) {
      throw std::invalid_argument("fscr::dot - Matrix cols != vector size");
    }
    
    std::vector<complex_d> res(mat.rows, complex_d());
    for (size_t i = 0; i < res.size(); ++i) {
      for (size_t j = 0; j < mat.cols; ++j) {
        res[i] += vec[j] * mat.mem[i * mat.cols + j];
      }
    }

    return res;
  }

  template<typename S>
  matrix_t<S> exp(const matrix_t<S>& mat) {
    auto res(mat);
    for (auto& m: res.mem) {
      m = std::exp(m);
    }
    return res;
  }

  template<typename S>
  std::vector<complex_d> expFlatten(const matrix_t<std::complex<S>> &mat) {
    std::vector<complex_d> res(mat.mem);
    for (auto& r: res) {
      r = std::exp(r);
    }
    return res;
  }

  /**
   * Section - Matrix type shorthands
   */
  using matrix = matrix_t<double>;
  using matrix_d = matrix_t<double>;
  using matrix_f = matrix_t<float>;
  using matrix_cd = matrix_t<complex_d>;
  using matrix_cf = matrix_t<std::complex<float>>;
}

#endif