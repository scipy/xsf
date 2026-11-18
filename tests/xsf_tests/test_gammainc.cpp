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

/*
If you want to check the expected values below with mpmath, this is the script:

import mpmath

mpmath.mp.dps = 50

for a, x in [(1e-300, 0.5), (1e-310, 0.5), (1e-300, 5.0)]:
    p = mpmath.gammainc(a, 0, x, regularized=True)
    q = mpmath.gammainc(a, x, mpmath.inf, regularized=True)
    print(f"({a!r}, {x!r}): gammainc = {float(p)!r}, gammaincc = {float(q)!r}")
*/

TEST_CASE("gammainc tiny a", "[gammainc][xsf_tests]") {
    // Regression test for gh-456: this exceeded 1 by 2e-14, or returned 0 for a < 1e-308.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [a, x, desired] = GENERATE(test_case{1e-300, 0.5, 1.0}, test_case{1e-310, 0.5, 1.0});

    const double out = xsf::gammainc(a, x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 2e-15;
    CAPTURE(a, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("gammaincc tiny a", "[gammaincc][xsf_tests]") {
    // Regression test for gh-456: relative error 2e-14.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [a, x, desired] =
        GENERATE(test_case{1e-300, 0.5, 5.597735947761608e-301}, test_case{1e-300, 5.0, 1.1482955912753257e-303});

    const double out = xsf::gammaincc(a, x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 2e-15;
    CAPTURE(a, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
