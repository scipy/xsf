#include "../testing_utils.h"

#include <xsf/cpu/bessel.h>
#include <xsf/mdspan.h>

TEST_CASE("rcty overflow", "[rcty][real][xsf_tests]") {
    // Regression test: x < 1e-60 returned -1e300 and 1e300, and the orders above |ry| > 1e300 were left unset.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 30

    # x * y_k(x)
    ry = lambda k, x: mpmath.sqrt(mpmath.pi * x / 2) * mpmath.bessely(k + mpmath.mpf(1) / 2, x)

    for k, x in [(2, "1e-70"), (149, 1), (1848, 1000), (200, 0.5)]:
        X = mpmath.mpf(x)
        dy = ry(k - 1, X) - k * ry(k, X) / X  # DLMF 10.51.2
        print(f"k={k}, x={x}: {float(ry(k, X))!r} {float(dy)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, int, double, double, double>;
    auto [n, k, x, desired_ry, desired_dy] = GENERATE(
        test_case{4, 2, 1e-70, -3e+140, 6e+210},
        test_case{150, 149, 1.0, -1.2573906731395004e+304, 1.873469766108334e+306},
        test_case{1848, 1848, 1000.0, -8.119236242111399e+307, 1.2616870528216194e+308},
        test_case{200, 200, 0.5, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}
    );

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> ry(n + 1, nan), dy(n + 1, nan);
    xsf::rcty(x, xsf::cxx::mdspan(ry.data(), ry.size()), xsf::cxx::mdspan(dy.data(), dy.size()));

    const double error_ry = xsf::extended_relative_error(ry[k], desired_ry);
    const double error_dy = xsf::extended_relative_error(dy[k], desired_dy);
    const double rtol = 1e-12;
    CAPTURE(n, k, x, ry[k], desired_ry, error_ry, dy[k], desired_dy, error_dy, rtol);
    REQUIRE(error_ry <= rtol);
    REQUIRE(error_dy <= rtol);
}

TEST_CASE("rcty n = 0", "[rcty][real][xsf_tests]") {
    // Regression test: with n = 0, ry[1] and dy[1] were written out of bounds.
    const double x = GENERATE(1e-70, 1.0);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> ry(2, nan), dy(2, nan);
    xsf::rcty(x, xsf::cxx::mdspan(ry.data(), 1), xsf::cxx::mdspan(dy.data(), 1));

    CAPTURE(x, ry[0], dy[0], ry[1], dy[1]);
    REQUIRE(ry[0] == -std::cos(x));
    REQUIRE(dy[0] == std::sin(x));
    REQUIRE(std::isnan(ry[1]));
    REQUIRE(std::isnan(dy[1]));
}
