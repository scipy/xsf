#include "../testing_utils.h"

#include <xsf/cephes/unity.h>

TEST_CASE("cephes::lgam1p near abs(x) = 1/2", "[cephes::lgam1p][xsf_tests]") {
    // Regression test for gh-348: the Taylor series stopped after 41 terms.
    // Reference values: float(mpmath.loggamma(1 + mpmath.mpf(x))) with mpmath.mp.dps = 50.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(
        test_case{-0.5, 0.5723649429247001}, test_case{0.5, -0.12078223763524522}, test_case{1.499, 0.2839799590505287}
    );

    const double out = xsf::cephes::lgam1p(x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-15;
    CAPTURE(x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
