cmake_minimum_required(VERSION 3.14)
project(xsf)

add_library(xsf INTERFACE)

# Scipy is restricted to C++17
# https://docs.scipy.org/doc/scipy/dev/toolchain.html#c-language-standards
target_compile_features(xsf INTERFACE cxx_std_17)

target_include_directories(xsf
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

install(
    TARGETS xsf
    EXPORT xsfTargets
)

install(
    DIRECTORY include/
    DESTINATION include
)

install(
    EXPORT xsfTargets
    FILE xsfTargets.cmake
    NAMESPACE xsf::
    DESTINATION lib/cmake/xsf
)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}" CACHE STRING "")
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}" CACHE STRING "")
endif()

# Tests
option(BUILD_TESTS "Build the tests" OFF)

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
