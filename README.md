# cpp-fft-fscr | C++ Fast Fourier Transform Library - From Scratch

**cpp-fft-fscr** is a header-only C++ library of [Discrete Fourier Transform (DFT)](https://en.wikipedia.org/wiki/Discrete_Fourier_transform) and  [Fast Fourier Transform (FFT)](https://en.wikipedia.org/wiki/Fast_Fourier_transform) functions created without any dependencies. This library is created From Scratch (fscr), literally.

In addition to provide a C++ library without dependencies, the **fscr** libraries are my attemp to make an easy-to-learn libraries for people who want to learn, especially, DFT and FFT implementation from scratch.

Features:
- A header-only library of DFT and FFT functions with various algorithms
- Matrix type created from scratch to support matrix operation in fft algorithms internal
- Written in C++11 format, and is C++11/14/17 compatible
- Released under a permissive, MIT license

To do:
- Prime-factor FFT algorithm - [wiki](https://en.wikipedia.org/wiki/Prime-factor_FFT_algorithm)
- Bruun's FFT algorithm - [wiki](https://en.wikipedia.org/wiki/Bruun%27s_FFT_algorithm)
- Rader's FFT algorithm - [wiki](https://en.wikipedia.org/wiki/Rader%27s_FFT_algorithm)
- Bluestein's FFT algorithm - [wiki](https://en.wikipedia.org/wiki/Bluestein%27s_FFT_algorithm)
- Hexagonal FFT - [wiki](https://en.wikipedia.org/wiki/Hexagonal_fast_Fourier_transform)
- SIMD technique for matrix operations
- More unit tests

## Algorithms
DFT:
- Naive looping
- Matrix operation

FFT:
- Divide and Conquer (Cooley-Tukey algorithm) - [wiki](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)

## Installation
Simply add the repository as the `git submodule` and/or add the header files to your project using:
``` C++
#include "cpp-fft-fscr/fft-fscr.hpp"
```

> Don't forget to add the installation directory into your include path

### CMake
Or, You can install the library from source using CMake.
``` bash
# clone cpp-fft-fscr from GitHub
git clone https://github.com/ikraduya/cpp-fft-fscr

# create build folder and specify installation location
cd cpp-fft-fscr
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/cpp-fft-fscr/installation/location

# install
cmake --install build

# Build tests (optional)
cmake -S . -B build_tests -DBUILD_TESTS=1
cmake --build build_tests --target tests
cd build_tests && ctest
```

For example, `/cpp-fft-fscr/installation/location` could be `/usr/local`.

If `DCMAKE_INSTALL_PREFIX` is not provided, the default installation location for linux would be `/usr/local`.

### Example Usage with CMake
You can see at `example_project/` directory on how to use it in your project. 

## Author
Ikraduya Edian

## License
MIT
