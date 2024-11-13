cmake_minimum_required(VERSION 3.14)
project(xsf)

# Scipy requ C++17
# https://docs.scipy.org/doc/scipy/dev/toolchain.html#c-language-standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Tests
enable_testing()
add_subdirectory(tests)
