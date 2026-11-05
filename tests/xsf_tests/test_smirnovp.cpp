#include "../testing_utils.h"

#include <xsf/cpu/stats.h>

// Regression: smirnovp returned -0.0 once n * (v - 1) overflowed int (n > 46340).
TEST_CASE("smirnovp large n", "[smirnovp][xsf_tests]") {
    const int n = GENERATE(46000, 50000, 100000, 1000000);
    const double x = 0.002;
    const double h = 1e-8;

    // Compare against a central difference of smirnov, which has no derivative terms.
    const double fd = (xsf::cpu::smirnov(n, x + h) - xsf::cpu::smirnov(n, x - h)) / (2 * h);
    const double out = xsf::cpu::smirnovp(n, x);
    CAPTURE(n, out, fd);
    REQUIRE(out < 0);
    REQUIRE(xsf::extended_relative_error(out, fd) <= 1e-6);
}
