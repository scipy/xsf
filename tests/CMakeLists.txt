if(NOT DEFINED ENV{XSREF_TABLES_PATH})
    message(FATAL_ERROR "XSREF_TABLES_PATH environment variable is not set. Specify the path to xsref/tables in your environment.")
endif()

set(XSREF_TABLES_PATH $ENV{XSREF_TABLES_PATH})

find_package(Catch2 3 REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)

add_library(xsf INTERFACE)
target_include_directories(xsf INTERFACE ${CMAKE_SOURCE_DIR}/include)

file(GLOB TEST_SOURCES "*/test_*.cpp")
foreach(test_file ${TEST_SOURCES})
  get_filename_component(test_name ${test_file} NAME_WE)
  add_executable(${test_name}.test ${test_file})
  target_link_libraries(${test_name}.test PRIVATE Catch2::Catch2WithMain Arrow::arrow_shared Parquet::parquet_shared xsf)
  target_compile_definitions(${test_name}.test PRIVATE XSREF_TABLES_PATH="${XSREF_TABLES_PATH}")
  include(CTest)
  include(Catch)
  catch_discover_tests(${test_name}.test)
endforeach()
