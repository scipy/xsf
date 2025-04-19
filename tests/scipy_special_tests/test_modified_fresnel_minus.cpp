#include "../testing_utils.h"

#include <xsf/fresnel.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "modified_fresnel_minus"};

TEST_CASE("modified_fresnel_minus d->DD scipy_special_tests", "[modified_fresnel_minus][d->DD][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<
            double, std::tuple<std::complex<double>, std::complex<double>, bool>, std::tuple<double, double>>(
            tables_path / "In_d-cd_cd.parquet", tables_path / "Out_d-cd_cd.parquet",
            tables_path / ("Err_d-cd_cd_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired0, desired1, fallback] = output;

    std::complex<double> out0;
    std::complex<double> out1;

    xsf::modified_fresnel_minus(x, out0, out1);
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
