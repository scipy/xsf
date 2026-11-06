#include "../testing_utils.h"

#include <xsf/cpu/struve.h>

TEST_CASE("it2struve0 moderate x", "[it2struve0][xsf_tests]") {
    // Regression test for gh-399 (https://github.com/scipy/xsf/issues/399):
    // the power series cancelled for moderate x, and the asymptotic expansion
    // was inaccurate for 24.5 <= x < 40.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 30

    for x in [24.45, 30.0]:
        # the integral of H0(t)/t from 0 to infinity is pi/2
        y = mpmath.pi / 2 - mpmath.quad(lambda t: mpmath.struveh(0, t) / t, mpmath.linspace(0, x, 31))
        print(f"it2struve0({x!r}) = {float(y)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(test_case{24.45, 0.026264940692299028}, test_case{30.0, 0.018152554625687087});

    const double out = xsf::it2struve0(x);
    const double error = xsf::extended_relative_error(out, desired);
    // pi/2 minus the integral from 0 to x cancels about 2 digits
    const double rtol = 2e-13;
    CAPTURE(x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
