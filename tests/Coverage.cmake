if(CMAKE_BUILD_TYPE STREQUAL "Coverage")

# Enable coverage compilation option
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /coverage")
endif()

# Add custom targets for generating coverage reports
add_custom_target(coverage
    COMMAND lcov --capture --directory . --output-file coverage.info
    COMMAND lcov --output-file coverage.info --extract coverage.info '*/include/xsf/*'
    COMMAND lcov --list coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating coverage report"
)

# Generate coverage reports in HTML format
add_custom_target(coverage_html
    COMMAND genhtml --demangle-cpp --legend coverage.info --output-directory coverage_report
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating HTML coverage report"
)
add_dependencies(coverage_html coverage)

endif() # CMAKE_BUILD_TYPE=Coverage