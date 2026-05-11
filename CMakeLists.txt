cmake_minimum_required(VERSION 3.14)
project(xsf)

# Scipy is restricted to C++17
# https://docs.scipy.org/doc/scipy/dev/toolchain.html#c-language-standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Tests
option(BUILD_TESTS "Build the tests" OFF)

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
