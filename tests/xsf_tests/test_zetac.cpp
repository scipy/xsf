#include "../testing_utils.h"

#include <xsf/zeta.h>

TEST_CASE("zetac large x", "[zetac][xsf_tests]") {
    // Regression test for gh-346: this returned 0 for 127 <= x < 1075.
    // Reference values: float(mpmath.zeta(x, 2)) with mpmath.mp.dps = 50.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(
        test_case{127.0, 5.877471754111438e-39}, test_case{1000.0, 9.332636185032189e-302}, test_case{1074.0, 5e-324}
    );

    const double out = xsf::zetac(x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 2e-16;
    CAPTURE(x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("zetac underflows to 0 for x >= 1076", "[zetac][xsf_tests]") {
    auto x = GENERATE(1076.0, 1e300);
    const double out = xsf::zetac(x);
    CAPTURE(x, out);
    REQUIRE(out == 0.0);
}
