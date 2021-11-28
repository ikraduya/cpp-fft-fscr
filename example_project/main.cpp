#include <vector>
#include <complex>
#include <iostream>
#include <limits>

#include "cpp-fft-fscr/fft-fscr.hpp"

template<typename T = double>
std::vector<T> generate_random_vec(size_t size, T lambda = 3.5)
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

template<typename T>
bool effectivelyEqual(T a, T b) {
  return std::abs(a - b) < std::numeric_limits<T>::epsilon();
}

template<typename T>
void printComplexVec(std::vector<std::complex<T>> vec) {
  std::cout << '[';
  for (const auto& v: vec) {
    std::cout << v << ',' << ' ';
  }
  std::cout << ']' << std::endl;
}

int main(int argc, char const *argv[]) {
  auto input_vector = generate_random_vec(1024);
  auto dft_loop_res = fscr::FFT::dft_loop(input_vector);
  auto dft_mat_res = fscr::FFT::dft_mat(input_vector);
  auto fft_dnc_res = fscr::FFT::fft_dnc(input_vector);

  auto less_complex_d = [](const std::complex<double>& lhs, const std::complex<double>& rhs){
    if (effectivelyEqual(lhs.real(), rhs.real())) {
      return lhs.imag() < rhs.imag();
    }
    return lhs.real() < rhs.real();
  };
  auto max_dft_loop = (*std::max_element(dft_loop_res.begin(), dft_loop_res.end(), less_complex_d));
  auto max_dft_mat = (*std::max_element(dft_mat_res.begin(), dft_mat_res.end(), less_complex_d));
  auto max_fft_dnc = (*std::max_element(fft_dnc_res.begin(), fft_dnc_res.end(), less_complex_d));

  std::cout << "Max from dft_loop = " << max_dft_loop << std::endl;
  std::cout << "Max from dft_mat = " << max_dft_mat << std::endl;
  std::cout << "Max from fft_dnc = " << max_fft_dnc << std::endl;

  return 0;
}