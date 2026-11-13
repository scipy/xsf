#include "../testing_utils.h"

#include <xsf/cephes/expn.h>

TEST_CASE("cephes::expn subnormal results", "[cephes::expn][xsf_tests]") {
    // Regression test: expn returned 0 for x > MAXLOG (about 709.78), but E_n(x)
    // is subnormal up to x of about 738.5 (n = 0).
    // Subnormal results are compared in units of the smallest subnormal.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    cases = [(50, 714.9428986597577), (1, 720.0), (0, 738.5281143251407),
             (52, 738.4605862096265), (69, 738.4394566275181)]
    for n, x in cases:
        print(f"expn({n}, {x!r}) = {float(mpmath.expint(n, x))!r}")
    */
    using test_case = std::tuple<int, double, double>;
    auto [n, x, desired] = GENERATE(
        test_case{50, 714.9428986597577, 4.174959247e-314}, test_case{1, 720.0, 2.81863344e-316},
        test_case{0, 738.5281143251407, 5e-324}, test_case{52, 738.4605862096265, 5e-324},
        test_case{69, 738.4394566275181, 0.0}
    );

    const double out = xsf::cephes::expn(n, x);
    // The exact values for n = 0, 52, 69 are just above, just above and just below half the smallest subnormal.
    const double tol = (desired < 1e-320) ? 0.0 : std::numeric_limits<double>::denorm_min();
    CAPTURE(n, x, out, desired, tol);
    REQUIRE(std::abs(out - desired) <= tol);
}
