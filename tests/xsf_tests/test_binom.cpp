#include "../testing_utils.h"

#include <xsf/binom.h>

TEST_CASE("binom k > 1e8 * abs(n)", "[binom][xsf_tests]") {
    // Regression test for gh-352: the asymptotic branch for large k missed
    // the n + 1 factor, overflowed an int cast for k >= 2**31, used sin
    // instead of sinpi, and returned 0 when k**n overflowed.

    /*
    import mpmath

    mpmath.mp.dps = 50

    for n, k in [(20.5, 2_100_000_000.25), (-10.5, 500_000_000_000.25), (1.0, 2_147_483_647.0),
                 (31.68727107709636, 10283293173.796093)]:
        print(f"binom({n!r}, {k!r}) = {float(mpmath.binomial(n, k))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [n, k, desired] = GENERATE(
        test_case{20.5, 2100000000.25, -9.316864069123419e-183},
        test_case{-10.5, 500000000000.25, 8.617145704224517e+104}, test_case{1.0, 2147483647.0, 0.0},
        test_case{31.68727107709636, 10283293173.796093, 5.0927883742793725e-294}
    );

    const double out = xsf::binom(n, k);
    const double error = xsf::extended_relative_error(out, desired);
    // limited by the two-term expansion and the rounding of k - floor(k) - n
    const double rtol = 2e-14;
    CAPTURE(n, k, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
