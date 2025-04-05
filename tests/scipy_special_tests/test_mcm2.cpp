#include "../testing_utils.h"

#include <xsf/mathieu.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "mcm2"};

TEST_CASE("mcm2 ddd->dd scipy_special_tests", "[mcm2][ddd->dd][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<
            std::tuple<double, double, double>, std::tuple<double, double, bool>, std::tuple<double, double>>(
            tables_path / "In_d_d_d-d_d.parquet", tables_path / "Out_d_d_d-d_d.parquet",
            tables_path / ("Err_d_d_d-d_d_" + get_platform_str() + ".parquet")
        )
    );

    auto [m, q, x] = input;
    auto [desired0, desired1, fallback] = output;

    double out0;
    double out1;

    xsf::mcm2(m, q, x, out0, out1);
    auto [tol0, tol1] = tol;

    auto error0 = xsf::extended_relative_error(out0, desired0);
    tol0 = adjust_tolerance(tol0);
    CAPTURE(m, q, x, out0, desired0, error0, tol0, fallback);
    REQUIRE(error0 <= tol0);

    auto error1 = xsf::extended_relative_error(out1, desired1);
    tol1 = adjust_tolerance(tol1);
    CAPTURE(m, q, x, out1, desired1, error1, tol1, fallback);
    REQUIRE(error1 <= tol1);
}
