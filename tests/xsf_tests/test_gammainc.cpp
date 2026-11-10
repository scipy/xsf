#include "../testing_utils.h"

#include <xsf/gamma.h>

/*
If you want to check the expected values below with mpmath, this is the script.
It prints log(x^a e^-x / Gamma(a)): this factor (and so gammainc if x < a, or
gammaincc if x > a) underflows to 0.

import mpmath

mpmath.mp.dps = 50

for a, x in [(1.7976931348623157e308, 3.0), (4.4e306, 3.7e180), (2.54e305, 1e308)]:
    print(f"({a!r}, {x!r}): {mpmath.mpf(a) * mpmath.log(x) - x - mpmath.loggamma(a)}")
*/

TEST_CASE("gammainc huge a", "[gammainc][gammaincc][xsf_tests]") {
    // Regression test for gh-462: these were NaN.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [a, x, desired] = GENERATE(
        test_case{std::numeric_limits<double>::max(), 3.0, 0.0}, test_case{4.4e306, 3.7e180, 0.0},
        test_case{2.54e305, 1e308, 1.0}
    );

    const double p = xsf::gammainc(a, x);
    const double q = xsf::gammaincc(a, x);
    CAPTURE(a, x, p, q, desired);
    REQUIRE(p == desired);
    REQUIRE(q == 1.0 - desired);
}
