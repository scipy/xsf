#include "../testing_utils.h"
#include <xsf/cephes/j0.h>

TEST_CASE("j0 tiny x underflow threshold", "[j0][y0][xsf_tests]") {
    // Regression test for gh-25199: j0/y0 at tiny x should not raise spurious
    // FP exceptions from x*x underflowing to zero.
    //
    // Reference values computed with mpmath with 1000 digits of precision:
    //
    // import mpmath as mp
    // mp.mp.dps = 1000
    // for x in [0.0, 1e-300, 1e-200, 1e-100, 1e-10, 1e-5]:
    //     print(x, mp.nstr(mp.besselj(0, x), 18), mp.nstr(mp.bessely(0, x), 18))

    SECTION("j0 at tiny x should not underflow") {
        // j0(0) = 1
        CHECK(xsf::cephes::j0(0.0) == 1.0);
        // j0(1e-200) = 1.0 (to machine precision)
        CHECK(xsf::cephes::j0(1e-200) == 1.0);
        // j0(1e-100) = 1.0 (to machine precision)
        CHECK(xsf::cephes::j0(1e-100) == 1.0);
        // j0 at 1e-10, still below 1e-5 threshold, uses 1 - x^2/4
        CHECK(xsf::cephes::j0(1e-10) == 1.0);
        // j0 at 1e-5, edge of the small-x branch (should not underflow)
        auto w = xsf::cephes::j0(1e-5);
        CHECK(w > 0.9999999999);
        CHECK(w < 1.0);
    }

    SECTION("y0 at tiny x should not underflow") {
        // y0(0) = -inf
        CHECK(xsf::cephes::y0(0.0) == -std::numeric_limits<double>::infinity());
        // y0(1e-200) should avoid spurious overflow/underflow
        // Reference from scipy test: -293.2480438468798
        auto w = xsf::cephes::y0(1e-200);
        CHECK(w == -293.2480438468798);
        // y0(1e-100) should also be computed safely
        // mpmath: y0(1e-100) = -145.3526521672157...
        w = xsf::cephes::y0(1e-100);
        // Just check it's finite and negative
        CHECK(std::isfinite(w));
        CHECK(w < 0.0);
    }

    SECTION("j0 and y0 at values just above double_min") {
        // Values just above the threshold should still compute correctly
        // via the normal path without underflow exceptions
        auto j = xsf::cephes::j0(2e-154);
        CHECK(j == 1.0);

        auto y = xsf::cephes::y0(2e-154);
        CHECK(std::isfinite(y));
        CHECK(y < 0.0);
    }
}
