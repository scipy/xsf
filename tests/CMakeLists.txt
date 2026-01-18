set(XSREF_TABLES_PATH "${CMAKE_SOURCE_DIR}/xsref/tables")
find_package(Catch2 3 REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)

#============================================
# SDB says: Use pkg-config to find blas & lapack
# which are requred to run the Mathieu tests.
# Find pkg-config
find_package(PkgConfig REQUIRED)

# Use pkg-config to find a library
pkg_search_module(BLAS REQUIRED IMPORTED_TARGET openblas blas)
pkg_search_module(LAPACK REQUIRED IMPORTED_TARGET openblas lapack)
#============================================

add_library(xsf INTERFACE)
target_include_directories(xsf INTERFACE ${CMAKE_SOURCE_DIR}/include)

set(TEST_BASE_DIR "${CMAKE_SOURCE_DIR}/tests")
include(${CMAKE_SOURCE_DIR}/tests/Coverage.cmake)

# Include CTest once at the top level
include(CTest)
enable_testing()

file(GLOB TEST_SOURCES "*/test_*.cpp")

foreach(test_file ${TEST_SOURCES})
  # Families of tests go in subfolders of xsf/tests. Test files in different
  # folders can have the same name. Try to generate a unique target name based
  # on the test name and its parent folder(s).
  get_filename_component(test_name ${test_file} NAME_WE)
  get_filename_component(test_dir ${test_file} DIRECTORY)
  file(RELATIVE_PATH test_dir ${TEST_BASE_DIR} ${test_dir})
  string(REPLACE "/" "-" test_dir ${test_dir})
  set(target_name ${test_dir}_${test_name})

  add_executable(${target_name} ${test_file})

  #===========================================
  # SDB added blas and lapack to link libs
  target_link_libraries(${target_name} PRIVATE
  Catch2::Catch2WithMain
  Arrow::arrow_shared
  Parquet::parquet_shared
  PkgConfig::BLAS
  PkgConfig::LAPACK
  xsf)

  target_compile_definitions(${target_name} PRIVATE XSREF_TABLES_PATH="${XSREF_TABLES_PATH}")
  include(CTest)

  target_compile_definitions(${target_name} PRIVATE 
    XSREF_TABLES_PATH="${XSREF_TABLES_PATH}"
  )
  
  # Manual test registration instead of catch_discover_tests
  add_test(
    NAME ${target_name}
    COMMAND ${target_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

  # Set test properties
  set_tests_properties(${target_name} PROPERTIES
    TIMEOUT 300
    ENVIRONMENT "XSREF_TABLES_PATH=${XSREF_TABLES_PATH}"
  )

endforeach()
