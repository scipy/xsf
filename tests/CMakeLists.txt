if(NOT DEFINED ENV{XSREF_TABLES_PATH})
    message(FATAL_ERROR "XSREF_TABLES_PATH environment variable is not set. Specify the path to xsref/tables in your environment.")
endif()

set(XSREF_TABLES_PATH $ENV{XSREF_TABLES_PATH})

find_package(Catch2 3 REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)

add_library(xsf INTERFACE)
target_include_directories(xsf INTERFACE ${CMAKE_SOURCE_DIR}/include)

set(TEST_BASE_DIR "${CMAKE_SOURCE_DIR}/tests")

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

  target_link_libraries(${target_name} PRIVATE Catch2::Catch2WithMain Arrow::arrow_shared Parquet::parquet_shared xsf)

  target_compile_definitions(${target_name} PRIVATE XSREF_TABLES_PATH="${XSREF_TABLES_PATH}")
  include(CTest)
  include(Catch)
  catch_discover_tests(${target_name})
endforeach()
