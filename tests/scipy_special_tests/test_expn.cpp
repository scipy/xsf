#include "../testing_utils.h"

#include <xsf/cephes/expn.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "expn"};

TEST_CASE("cephes::expn dd->d scipy_special_tests", "[cephes::expn][dd->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(xsf_test_cases<std::tuple<double, double>, std::tuple<double, bool>, double>(
        tables_path / "In_d_d-d.parquet", tables_path / "Out_d_d-d.parquet",
        tables_path / ("Err_d_d-d_" + get_platform_str() + ".parquet")
    ));

    auto [n, x] = input;
    auto [desired, fallback] = output;
    auto out = xsf::cephes::expn(n, x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(n, x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}

TEST_CASE("cephes::expn pd->d scipy_special_tests", "[cephes::expn][pd->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] =
        GENERATE(xsf_test_cases<std::tuple<std::ptrdiff_t, double>, std::tuple<double, bool>, double>(
            tables_path / "In_p_d-d.parquet", tables_path / "Out_p_d-d.parquet",
            tables_path / ("Err_p_d-d_" + get_platform_str() + ".parquet")
        ));

    auto [n, x] = input;
    auto [desired, fallback] = output;
    auto out = xsf::cephes::expn(n, x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(n, x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
