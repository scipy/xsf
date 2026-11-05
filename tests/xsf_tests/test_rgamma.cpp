#include "../testing_utils.h"

#include <xsf/cephes/rgamma.h>

TEST_CASE("cephes::rgamma -178 < x < -170", "[cephes::rgamma][xsf_tests]") {
    // Regression test for gh-377: this returned inf because Gamma(x) underflows.
    // Reference values: float(mpmath.rgamma(x)) with mpmath.mp.dps = 50.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(
        test_case{-170.91716291792244, -6.637711907120738e+307}, test_case{-175.99999999999997, 5.624744700044757e+306},
        test_case{-176.00000000000003, -5.624744700046411e+306}
    );

    const double out = xsf::cephes::rgamma(x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-15;
    CAPTURE(x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("cephes::rgamma tiny negative subnormal x", "[cephes::rgamma][xsf_tests]") {
    // Regression test for gh-417: this returned -0.0 because 1/x overflows.
    // Reference values: float(mpmath.rgamma(x)) with mpmath.mp.dps = 50.
    SET_FP_FORMAT()
    REQUIRE(xsf::cephes::rgamma(-1e-310) == -1e-310);
    REQUIRE(xsf::cephes::rgamma(-5e-324) == -5e-324);
}
