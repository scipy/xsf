#include "../testing_utils.h"

#include <xsf/erf.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "erf"};

TEST_CASE("erf f->f scipy_special_tests", "[erf][f->f][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<float, std::tuple<float, bool>, float>(
            tables_path / "In_f-f.parquet", tables_path / "Out_f-f.parquet",
            tables_path / ("Err_f-f_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired, fallback] = output;
    auto out = xsf::erf(x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}

TEST_CASE("erf D->D scipy_special_tests", "[erf][D->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<std::complex<double>, std::tuple<std::complex<double>, bool>, double>(
            tables_path / "In_cd-cd.parquet", tables_path / "Out_cd-cd.parquet",
            tables_path / ("Err_cd-cd_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired, fallback] = output;
    auto out = xsf::erf(x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}

TEST_CASE("erf d->d scipy_special_tests", "[erf][d->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<double, std::tuple<double, bool>, double>(
            tables_path / "In_d-d.parquet", tables_path / "Out_d-d.parquet",
            tables_path / ("Err_d-d_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired, fallback] = output;
    auto out = xsf::erf(x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
