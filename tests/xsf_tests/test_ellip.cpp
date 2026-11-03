#include "../testing_utils.h"

#include <xsf/ellip.h>

TEST_CASE("ellipkinc and ellipeinc huge phi", "[ellipkinc][ellipeinc][xsf_tests]") {
    // Regression test for gh-359: both returned NaN.
    // Reference values: float(mpmath.ellipf(phi, m)), float(mpmath.ellipe(phi, m)) with mpmath.mp.dps = 50.
    SET_FP_FORMAT()
    const double phi = -1.7885758944711641e19, m = -1.2443495040003017e35;
    const double rtol = 2e-16;

    const double outk = xsf::ellipkinc(phi, m);
    const double errk = xsf::extended_relative_error(outk, -1348.954827564235);
    CAPTURE(outk, errk, rtol);
    REQUIRE(errk <= rtol);

    const double oute = xsf::ellipeinc(phi, m);
    const double erre = xsf::extended_relative_error(oute, -4.0166009415816355e+36);
    CAPTURE(oute, erre, rtol);
    REQUIRE(erre <= rtol);

    // n * E(m) overflows here, so this must stay -inf, not NaN
    REQUIRE(xsf::ellipeinc(-7.49757919668992e261, -2.6578197111983853e169) == -std::numeric_limits<double>::infinity());
}
