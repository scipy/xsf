#include "../testing_utils.h"

#include <xsf/binom.h>

// Regression: for k > 1e8 |n|, binom took dk = k - floor(k) and the parity of floor(k) only when
// floor(k) fit in an int. Beyond INT_MAX it used sin(k pi) directly, which loses all accuracy for
// large k, so e.g. binom(0.5, 3000000000.5) was -1.8e-21 instead of 0.
TEST_CASE("binom large k", "[binom][xsf_tests]") {
    // 1 + n - k is a non-positive integer here, so the result is exactly 0.
    REQUIRE(xsf::binom(0.5, 3000000000.5) == 0.0);

    // Reference values from mpmath.binomial.
    const double cases[][3] = {
        {0.5, 3e9, -1.7167742314287727e-15},
        {0.5, 1e12 + 1, 2.8209479177356079e-19},
        {0.5, 1e15, -8.9206205807638589e-24},
        {-0.5, 1e15, 1.7841241161527709e-8},
    };
    for (const auto &c : cases) {
        const double out = xsf::binom(c[0], c[1]);
        CAPTURE(c[0], c[1], out, c[2]);
        REQUIRE(xsf::extended_relative_error(out, c[2]) <= 1e-6);
    }
}
