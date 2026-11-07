#include "../testing_utils.h"

#include <xsf/cpu/fresnel.h>

TEST_CASE("fresnel huge x", "[fresnel][xsf_tests]") {
    // Regression test for gh-386: this returned NaN for |x| >= 2^512, where x * x overflows.
    // mpmath.fresnels(x) and mpmath.fresnelc(x) round to +-0.5 here.
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(test_case{0x1p512, 0.5}, test_case{1e155, 0.5}, test_case{-1e200, -0.5});

    double s, c;
    xsf::fresnel(x, s, c);
    CAPTURE(x, s, c, desired);
    REQUIRE(s == desired);
    REQUIRE(c == desired);
}
