if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  set(COVERAGE_DIR share/xsf/cov)

  # Enable coverage compilation option
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
  endif()

  # Number of path components in the build dir, needed for GCOV_PREFIX_STRIP
  # at test time so .gcda files land in share/xsf/cov/ mirroring the .gcno layout
  string(REGEX MATCHALL "/" _slashes "${CMAKE_BINARY_DIR}")
  list(LENGTH _slashes GCOV_PREFIX_STRIP_COUNT)

  # Install .gcno files (compile-time output) into the package, preserving
  # their relative layout so they line up with the .gcda files written later
  install(DIRECTORY ${CMAKE_BINARY_DIR}/
    DESTINATION ${COVERAGE_DIR}/
    FILES_MATCHING PATTERN "*.gcno"
  )

  # Drop a small file recording the strip count
  file(WRITE ${CMAKE_BINARY_DIR}/gcov_prefix_strip.txt "${GCOV_PREFIX_STRIP_COUNT}")
  install(FILES ${CMAKE_BINARY_DIR}/gcov_prefix_strip.txt DESTINATION ${COVERAGE_DIR}/)
endif() # CMAKE_BUILD_TYPE=Coverage
