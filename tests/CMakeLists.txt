find_package(Catch2 3 REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)

add_library(xsf INTERFACE)

target_include_directories(xsf INTERFACE ${CMAKE_SOURCE_DIR}/include)

file(GLOB TEST_SOURCES "test_*.cpp")

foreach(test_file in LISTS TEST_SOURCES)
  get_filename_component(test_name ${test_file} NAME_WE)
  add_executable(${test_name} ${test_file})
  target_link_libraries(${test_name} PRIVATE Catch2::Catch2WithMain Arrow::arrow_shared Parquet::parquet_shared xsf)

  include(CTest)
  include(Catch)
  catch_discover_tests(${test_name})
endforeach()
