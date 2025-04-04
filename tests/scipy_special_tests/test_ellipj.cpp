#include "../testing_utils.h"

#include <xsf/ellip.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "ellipj"};

TEST_CASE("ellipj dd->dddd scipy_special_tests", "[ellipj][dd->dddd][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(xsf_test_cases<
                                         std::tuple<double, double>, std::tuple<double, double, double, double, bool>,
                                         std::tuple<double, double, double, double>>(
        tables_path / "In_d_d-d_d_d_d.parquet", tables_path / "Out_d_d-d_d_d_d.parquet",
        tables_path / ("Err_d_d-d_d_d_d_" + get_platform_str() + ".parquet")
    ));

    auto [u, m] = input;
    auto [desired0, desired1, desired2, desired3, fallback] = output;

    double out0;
    double out1;
    double out2;
    double out3;

    xsf::ellipj(u, m, out0, out1, out2, out3);
    auto [tol0, tol1, tol2, tol3] = tol;

    auto error0 = xsf::extended_relative_error(out0, desired0);
    tol0 = adjust_tolerance(tol0);
    CAPTURE(u, m, out0, desired0, error0, tol0, fallback);
    REQUIRE(error0 <= tol0);

    auto error1 = xsf::extended_relative_error(out1, desired1);
    tol1 = adjust_tolerance(tol1);
    CAPTURE(u, m, out1, desired1, error1, tol1, fallback);
    REQUIRE(error1 <= tol1);

    auto error2 = xsf::extended_relative_error(out2, desired2);
    tol2 = adjust_tolerance(tol2);
    CAPTURE(u, m, out2, desired2, error2, tol2, fallback);
    REQUIRE(error2 <= tol2);

    auto error3 = xsf::extended_relative_error(out3, desired3);
    tol3 = adjust_tolerance(tol3);
    CAPTURE(u, m, out3, desired3, error3, tol3, fallback);
    REQUIRE(error3 <= tol3);
}
