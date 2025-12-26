cmake_minimum_required(VERSION 3.14)
project(xsf)

# Scipy is restricted to C++17
# https://docs.scipy.org/doc/scipy/dev/toolchain.html#c-language-standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Help find packages in conda/pixi environment
if(DEFINED ENV{CONDA_PREFIX})
    message(STATUS "Using conda prefix: $ENV{CONDA_PREFIX}")
    list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
    set(ENV{PKG_CONFIG_PATH} "$ENV{CONDA_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()

#============================================
# SDB says: Use pkg-config to find blas & lapack
# which are requred to run the Mathieu tests.
# Find pkg-config
find_package(PkgConfig REQUIRED)

# Use pkg-config to find a library
#pkg_check_modules(BLAS REQUIRED blas)
#pkg_check_modules(LAPACK REQUIRED lapack)

#pkg_search_module(BLAS IMPORTED_TARGET REQUIRED openblas blas)
#pkg_search_module(LAPACK IMPORTED_TARGET REQUIRED openblas lapack)

#============================================

# Tests
# Turned on by SDB
option(BUILD_TESTS "Build the tests" ON) 

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
