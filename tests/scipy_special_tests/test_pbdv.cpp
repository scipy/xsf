#include "../testing_utils.h"

#include <xsf/par_cyl.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "pbdv"};

TEST_CASE("pbdv dd->dd scipy_special_tests", "[pbdv][dd->dd][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<std::tuple<double, double>, std::tuple<double, double, bool>, std::tuple<double, double>>(
            tables_path / "In_d_d-d_d.parquet", tables_path / "Out_d_d-d_d.parquet",
            tables_path / ("Err_d_d-d_d_" + get_platform_str() + ".parquet")
        )
    );

    auto [v, x] = input;
    auto [desired0, desired1, fallback] = output;

    double out0;
    double out1;

    xsf::pbdv(v, x, out0, out1);
    auto [tol0, tol1] = tol;

    auto error0 = xsf::extended_relative_error(out0, desired0);
    tol0 = adjust_tolerance(tol0);
    CAPTURE(v, x, out0, desired0, error0, tol0, fallback);
    REQUIRE(error0 <= tol0);

    auto error1 = xsf::extended_relative_error(out1, desired1);
    tol1 = adjust_tolerance(tol1);
    CAPTURE(v, x, out1, desired1, error1, tol1, fallback);
    REQUIRE(error1 <= tol1);
}

TEST_CASE("pbdv remains stable for large arguments", "[pbdv][regression]") {
    struct Case {
        double v;
        double x;
        double desired0;
        double desired1;
    };
    const Case cases[] = {
        {std::nextafter(50.0, 0.0), -std::sqrt(500.0), 3.1429284557854132e37, -2.7008764338749622e38},
        {std::nextafter(50.0, 0.0), -std::nextafter(50.0, 0.0), 5.6977287669133012e235, -1.3650529001419031e237},
        {-1.0, 50.0, 7.3587705514938545e-274, -1.8411632169283370e-272},
    };

    for (const auto &test_case : cases) {
        double out0;
        double out1;
        xsf::pbdv(test_case.v, test_case.x, out0, out1);
        CAPTURE(test_case.v, test_case.x, out0, out1);
        REQUIRE(xsf::extended_relative_error(out0, test_case.desired0) <= 5e-13);
        REQUIRE(xsf::extended_relative_error(out1, test_case.desired1) <= 5e-13);
    }
}
