#include "../testing_utils.h"

#include <xsf/cephes/chdtr.h>

TEST_CASE("cephes::chdtr subnormal x", "[cephes::chdtr][xsf_tests]") {
    // Regression test for gh-541 (https://github.com/scipy/xsf/issues/541):
    // x / 2 is subnormal (or 0) for x < 2 * DBL_MIN, which lost precision,
    // e.g. chdtr(1, 5e-324) was 0.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for df, x in [(1, 5e-324), (1, 1.5e-323), (1, 1e-310), (0.5, 1e-310), (2.002, 3 * 5e-324)]:
        v = mpmath.gammainc(mpmath.mpf(df) / 2, 0, mpmath.mpf(x) / 2, regularized=True)
        print(f"chdtr({df!r}, {x!r}) = {float(v)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [df, x, desired] = GENERATE(
        test_case{1, 5e-324, 1.7735048886036274e-162}, test_case{1, 1.5e-323, 3.0718005745332645e-162},
        test_case{1, 1e-310, 7.978845608028641e-156}, test_case{0.5, 1e-310, 2.933738615886127e-78},
        test_case{2.002, 3 * 5e-324, 5e-324}
    );

    const double out = xsf::cephes::chdtr(df, x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-15;
    CAPTURE(df, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("cephes::chdtr cephes::chdtrc x = 0 and subnormal df", "[cephes::chdtr][cephes::chdtrc][xsf_tests]") {
    // Regression test for gh-541: df / 2 underflows to 0, and igam(0, 0) and igamc(0, 0) are NaN,
    // but the CDF is 0 at x = 0 for any df > 0 (DLMF 8.2.4).
    const double df = 5e-324;
    REQUIRE(xsf::cephes::chdtr(df, 0.0) == 0.0);
    REQUIRE(xsf::cephes::chdtrc(df, 0.0) == 1.0);
}
