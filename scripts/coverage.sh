#!/usr/bin/env sh
export GCOV_PREFIX_STRIP=$(cat "$CONDA_PREFIX/share/xsf/cov/gcov_prefix_strip.txt")
lcov --capture --directory "$COVERAGE_DIR" --output-file "$COVERAGE_DIR/coverage.info" --exclude "*/catch2/*" --exclude "*/arrow/*" --exclude "*/parquet/*"
lcov --output-file "$COVERAGE_DIR/coverage.info" --extract "$COVERAGE_DIR/coverage.info" "*/include/xsf/*"
lcov --list "$COVERAGE_DIR/coverage.info"
genhtml --demangle-cpp --legend "$COVERAGE_DIR/coverage.info" --output-directory "$COVERAGE_DIR/coverage_report"
