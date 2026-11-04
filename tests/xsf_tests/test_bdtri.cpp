#include "../testing_utils.h"

#include <xsf/stats.h>

TEST_CASE("bdtri small y and k >= n / 2", "[bdtri][xsf_tests]") {
    // Regression test for gh-353: 1 - y lost the low digits of small y.

    /*
    Reference values from mpmath:

    import mpmath

    mpmath.mp.dps = 50

    for k, n, y in [(96, 191, 2.1e-22), (96, 191, 1e-12)]:
        # bdtr(k, n, p) = I_{1-p}(n - k, k + 1)
        f = lambda p: mpmath.log(mpmath.betainc(n - k, k + 1, 0, 1 - p, regularized=True) / y)
        p = mpmath.findroot(f, (0.5, 0.9), solver="anderson")
        print(f"bdtri({k}, {n}, {y!r}) = {float(p)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, int, double, double>;
    auto [k, n, y, desired] =
        GENERATE(test_case{96, 191, 2.1e-22, 0.8144065639961107}, test_case{96, 191, 1e-12, 0.7429714216811626});

    const double out = xsf::bdtri(k, n, y);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 2e-16;
    CAPTURE(k, n, y, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
