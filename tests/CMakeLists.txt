set(Python3_FIND_VIRTUALENV ONLY)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import xsref.tables as t; print(t.get_tables_path())"
  OUTPUT_VARIABLE PY_OUT
  ERROR_VARIABLE PY_ERR
  RESULT_VARIABLE PY_RES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT PY_RES EQUAL 0)
  message(FATAL_ERROR "Failed to get xsref tables path: ${PY_ERR}")
endif()

set(XSREF_TABLES_PATH "${PY_OUT}" CACHE PATH "Path to xsref tables" FORCE)

find_package(Catch2 3 REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # Enable debug assertions to make sure they won't get triggered in debug builds of
  # downstream code.
    if(MSVC OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_SIMULATE_ID MATCHES "MSVC"))
        add_compile_definitions(_ITERATOR_DEBUG_LEVEL=2)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        add_compile_definitions(_GLIBCXX_ASSERTIONS)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_definitions(_LIBCPP_ENABLE_ASSERTIONS=1)
    endif()
endif()

set(TEST_BASE_DIR "${CMAKE_SOURCE_DIR}/tests")

include(${CMAKE_SOURCE_DIR}/tests/Coverage.cmake)

file(GLOB_RECURSE TEST_SOURCES "${TEST_BASE_DIR}/*/test_*.cpp")
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
