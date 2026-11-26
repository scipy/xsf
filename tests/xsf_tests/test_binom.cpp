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

TEST_CASE("binom large non-integer abs(k)", "[binom][xsf_tests]") {
    // Regression test for gh-352: the beta fallback lost ~abs(k) eps for
    // large non-integer abs(k). The third case over- or underflows poch, and
    // the last one has k - n close to an integer.

    /*
    import mpmath

    mpmath.mp.dps = 50

    for n, k in [(31.69, 1_000_000_000.796), (-7.25, -1_000_000_000_000.5),
                 (-38.69913174872936, 2169303697.7401443), (0.49999999999999994, 20000.5)]:
        print(f"binom({n!r}, {k!r}) = {float(mpmath.binomial(n, k))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [n, k, desired] = GENERATE(
        test_case{31.69, 1000000000.796, -5.746138060441207e-261},
        test_case{-7.25, -1000000000000.5, -1.224023543149038e+72},
        test_case{-38.69913174872936, 2169303697.7401443, 6.492179234150414e+307},
        test_case{0.49999999999999994, 20000.5, 1.73929015396726e-23}
    );

    const double out = xsf::binom(n, k);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-14;
    CAPTURE(n, k, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
