find_package(Catch2 3 REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)


add_library(xsf INTERFACE)
target_include_directories(xsf INTERFACE ${CMAKE_SOURCE_DIR}/include)

add_executable(test_hyp2f1.o test_hyp2f1.cpp)

target_link_libraries(test_hyp2f1.o PRIVATE Catch2::Catch2WithMain Arrow::arrow_shared Parquet::parquet_shared xsf)
include(CTest)
include(Catch)
catch_discover_tests(test_hyp2f1.o)
