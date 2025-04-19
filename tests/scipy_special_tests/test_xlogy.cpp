#include "../testing_utils.h"

#include <xsf/log.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "xlogy"};

TEST_CASE("xlogy dd->d scipy_special_tests", "[xlogy][dd->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<std::tuple<double, double>, std::tuple<double, bool>, double>(
            tables_path / "In_d_d-d.parquet", tables_path / "Out_d_d-d.parquet",
            tables_path / ("Err_d_d-d_" + get_platform_str() + ".parquet")
        )
    );

    auto [x, y] = input;
    auto [desired, fallback] = output;
    auto out = xsf::xlogy(x, y);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, y, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}

TEST_CASE("xlogy DD->D scipy_special_tests", "[xlogy][DD->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<
            std::tuple<std::complex<double>, std::complex<double>>, std::tuple<std::complex<double>, bool>, double>(
            tables_path / "In_cd_cd-cd.parquet", tables_path / "Out_cd_cd-cd.parquet",
            tables_path / ("Err_cd_cd-cd_" + get_platform_str() + ".parquet")
        )
    );

    auto [x, y] = input;
    auto [desired, fallback] = output;
    auto out = xsf::xlogy(x, y);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, y, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}

TEST_CASE("xlogy ff->f scipy_special_tests", "[xlogy][ff->f][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<std::tuple<float, float>, std::tuple<float, bool>, float>(
            tables_path / "In_f_f-f.parquet", tables_path / "Out_f_f-f.parquet",
            tables_path / ("Err_f_f-f_" + get_platform_str() + ".parquet")
        )
    );

    auto [x, y] = input;
    auto [desired, fallback] = output;
    auto out = xsf::xlogy(x, y);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, y, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
