# cpp-fft-fscr Example Project

This is an example on how to use **cpp-fft-fscr** library.

Make sure that cpp-fft-fscr library is already installed.

## Usage
``` bash
cmake -S . -B build -DCMAKE_PREFIX_PATH='/cpp-fft-fscr/installation/location'
cmake --build build
```

In `CMakeLists.txt`, we just need to add
``` cmake
find_package(cpp-fft-fscr CONFIG REQUIRED)
if (cpp-fft-fscr_FOUND)
  message(STATUS "cpp-fft-fscr library is found")
endif()

target_link_libraries(main PRIVATE cpp-fft-fscr::cpp-fft-fscr) # installed cpp-fft-fscr include/ path is automatically added
``` 
