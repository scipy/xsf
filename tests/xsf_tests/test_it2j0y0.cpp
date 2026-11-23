#include "../testing_utils.h"

#include <xsf/cpu/bessel.h>

TEST_CASE("it2j0y0 moderate x", "[it2j0y0][xsf_tests]") {
    // Regression test for gh-446 (https://github.com/scipy/xsf/issues/446):
    // the power series cancelled for moderate x, and the asymptotic expansion
    // was inaccurate for 20 < x < 34.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 30

    L = lambda t: mpmath.log(t / 2) + mpmath.euler
    f = lambda t: (mpmath.bessely(0, t) - 2 / mpmath.pi * L(t)) / t
    for x in [19.5, 21.28]:
        ts = mpmath.linspace(0, x, 21)
        j0int = mpmath.quad(lambda t: (1 - mpmath.besselj(0, t)) / t, ts)
        # int_x^inf Y0(t)/t dt = pi/6 - L(x)^2/pi - int_0^x f(t) dt
        y0int = mpmath.pi / 6 - L(x) ** 2 / mpmath.pi - mpmath.quad(f, ts)
        print(f"it2j0y0({x!r}) = {float(j0int)!r}, {float(y0int)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, desired_j0, desired_y0] = GENERATE(
        test_case{19.5, 2.856464784421553, 0.008984103411895046},
        test_case{21.28, 2.933756950893377, 9.723256856027431e-06}
    );

    double j0int, y0int;
    xsf::it2j0y0(x, j0int, y0int);
    const double error_j0 = xsf::extended_relative_error(j0int, desired_j0);
    const double error_y0 = xsf::extended_relative_error(y0int, desired_y0);
    // the Y0 integral is small and cancels about 3 digits, and more near its zero at x = 21.28
    const double rtol_j0 = 1e-14, rtol_y0 = 1e-9;
    CAPTURE(x, j0int, desired_j0, error_j0, y0int, desired_y0, error_y0);
    REQUIRE(error_j0 <= rtol_j0);
    REQUIRE(error_y0 <= rtol_y0);
}
