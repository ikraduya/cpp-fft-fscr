cmake_minimum_required(VERSION 3.14)

project(Example_cpp-fft-fscr)

find_package(cpp-fft-fscr CONFIG REQUIRED)
if (cpp-fft-fscr_FOUND)
  message(STATUS "cpp-fft-fscr library is found")
elseif()
  message(FATAL_ERROR "cpp-fft-fscr library is not found")
endif()

# An executable for the file(s)
add_executable(main main.cpp)
target_link_libraries(main PRIVATE cpp-fft-fscr::cpp-fft-fscr)

# Set the properties required, e.g. C++ standard used. Apply to library and examples
set_target_properties(
  main
  PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
)