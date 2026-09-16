#include "../testing_utils.h"
#include <cfenv>
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
        using test_case = std::tuple<double, double>;
        auto [x, expected] = GENERATE(
            test_case{0.0, 1.0}, test_case{1e-200, 1.0}, test_case{1e-100, 1.0}, test_case{1e-10, 1.0},
            test_case{2e-154, 1.0}
        );
        std::feclearexcept(FE_ALL_EXCEPT);
        CHECK(xsf::cephes::j0(x) == expected);
        CHECK(std::fetestexcept(FE_UNDERFLOW) == 0);
    }

    SECTION("j0 at 1e-5 uses 1 - x^2/4") {
        std::feclearexcept(FE_ALL_EXCEPT);
        auto w = xsf::cephes::j0(1e-5);
        CHECK(w > 0.9999999999);
        CHECK(w < 1.0);
        CHECK(std::fetestexcept(FE_UNDERFLOW) == 0);
    }

    SECTION("y0 at tiny x should not underflow") {
        using test_case = std::tuple<double, double>;
        auto [x, expected] =
            GENERATE(test_case{0.0, -std::numeric_limits<double>::infinity()}, test_case{1e-200, -293.2480438468798});
        std::feclearexcept(FE_ALL_EXCEPT);
        CHECK(xsf::cephes::y0(x) == expected);
        CHECK(std::fetestexcept(FE_UNDERFLOW) == 0);
    }

    SECTION("y0 at tiny x is finite with no underflow") {
        auto x = GENERATE(1e-100, 2e-154);
        std::feclearexcept(FE_ALL_EXCEPT);
        auto w = xsf::cephes::y0(x);
        CHECK(std::isfinite(w));
        CHECK(w < 0.0);
        CHECK(std::fetestexcept(FE_UNDERFLOW) == 0);
    }
}
