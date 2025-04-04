#include "../testing_utils.h"

#include <xsf/kelvin.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "kelvin"};

TEST_CASE("kelvin d->DDDD scipy_special_tests", "[kelvin][d->DDDD][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<
            double,
            std::tuple<std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>, bool>,
            std::tuple<double, double, double, double>>(
            tables_path / "In_d-cd_cd_cd_cd.parquet", tables_path / "Out_d-cd_cd_cd_cd.parquet",
            tables_path / ("Err_d-cd_cd_cd_cd_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired0, desired1, desired2, desired3, fallback] = output;

    std::complex<double> out0;
    std::complex<double> out1;
    std::complex<double> out2;
    std::complex<double> out3;

    xsf::kelvin(x, out0, out1, out2, out3);
    auto [tol0, tol1, tol2, tol3] = tol;

    auto error0 = xsf::extended_relative_error(out0, desired0);
    tol0 = adjust_tolerance(tol0);
    CAPTURE(x, out0, desired0, error0, tol0, fallback);
    REQUIRE(error0 <= tol0);

    auto error1 = xsf::extended_relative_error(out1, desired1);
    tol1 = adjust_tolerance(tol1);
    CAPTURE(x, out1, desired1, error1, tol1, fallback);
    REQUIRE(error1 <= tol1);

    auto error2 = xsf::extended_relative_error(out2, desired2);
    tol2 = adjust_tolerance(tol2);
    CAPTURE(x, out2, desired2, error2, tol2, fallback);
    REQUIRE(error2 <= tol2);

    auto error3 = xsf::extended_relative_error(out3, desired3);
    tol3 = adjust_tolerance(tol3);
    CAPTURE(x, out3, desired3, error3, tol3, fallback);
    REQUIRE(error3 <= tol3);
}
