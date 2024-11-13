# Get googletest
# https://google.github.io/googletest/quickstart-cmake.html
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/releases/download/v1.15.2/googletest-1.15.2.tar.gz
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

# Add Tests
function(add_gtest test_name)
    add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name} GTest::gtest_main)
    gtest_discover_tests(${test_name})
endfunction()

add_gtest(airy)
