#include "../testing_utils.h"

#include <xsf/gamma.h>

TEST_CASE("gammaln subnormal arguments", "[gammaln][xsf_tests]") {
    // Regression test for gh-540: the x < 13 branch divided by x, overflowed, and returned inf.
    // Reference values: float(mpmath.log(abs(mpmath.gamma(x)))) with mpmath.mp.dps = 60.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(test_case{5e-309, 709.889355822726}, test_case{-5e-309, 709.889355822726});

    const double out = xsf::gammaln(x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-16;
    CAPTURE(x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
