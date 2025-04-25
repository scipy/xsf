#include "../testing_utils.h"

#include <xsf/sici.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "shichi"};

TEST_CASE("shichi d->dd scipy_special_tests", "[shichi][d->dd][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<double, std::tuple<double, double, bool>, std::tuple<double, double>>(
            tables_path / "In_d-d_d.parquet", tables_path / "Out_d-d_d.parquet",
            tables_path / ("Err_d-d_d_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired0, desired1, fallback] = output;

    double out0;
    double out1;

    xsf::shichi(x, out0, out1);
    auto [tol0, tol1] = tol;

    auto error0 = xsf::extended_relative_error(out0, desired0);
    tol0 = adjust_tolerance(tol0);
    CAPTURE(x, out0, desired0, error0, tol0, fallback);
    REQUIRE(error0 <= tol0);

    auto error1 = xsf::extended_relative_error(out1, desired1);
    tol1 = adjust_tolerance(tol1);
    CAPTURE(x, out1, desired1, error1, tol1, fallback);
    REQUIRE(error1 <= tol1);
}

TEST_CASE("shichi D->DD scipy_special_tests", "[shichi][D->DD][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<
            std::complex<double>, std::tuple<std::complex<double>, std::complex<double>, bool>,
            std::tuple<double, double>>(
            tables_path / "In_cd-cd_cd.parquet", tables_path / "Out_cd-cd_cd.parquet",
            tables_path / ("Err_cd-cd_cd_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired0, desired1, fallback] = output;

    std::complex<double> out0;
    std::complex<double> out1;

    xsf::shichi(x, out0, out1);
    auto [tol0, tol1] = tol;

    auto error0 = xsf::extended_relative_error(out0, desired0);
    tol0 = adjust_tolerance(tol0);
    CAPTURE(x, out0, desired0, error0, tol0, fallback);
    REQUIRE(error0 <= tol0);

    auto error1 = xsf::extended_relative_error(out1, desired1);
    tol1 = adjust_tolerance(tol1);
    CAPTURE(x, out1, desired1, error1, tol1, fallback);
    REQUIRE(error1 <= tol1);
}
