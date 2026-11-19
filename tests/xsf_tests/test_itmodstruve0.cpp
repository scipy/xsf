#include "../testing_utils.h"

#include <xsf/cpu/struve.h>

TEST_CASE("itmodstruve0 20 < x < 40 and x > 40", "[itmodstruve0][xsf_tests]") {
    // Regression test for gh-400 (https://github.com/scipy/xsf/issues/400):
    // the asymptotic expansion was used for too small x.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 30

    for x in [21.3, 30.0, 45.0]:
        y = mpmath.quad(lambda t: mpmath.struvel(0, t), [0, x])
        print(f"itmodstruve0({x!r}) = {float(y)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(
        test_case{21.3, 158791117.07105565}, test_case{30.0, 795538858181.7499}, test_case{45.0, 2.1075225881890097e+18}
    );

    const double out = xsf::itmodstruve0(x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-15;
    CAPTURE(x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
