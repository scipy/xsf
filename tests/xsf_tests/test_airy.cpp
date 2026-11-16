#include "../testing_utils.h"

#include <xsf/cpu/airy.h>

TEST_CASE("airy overflow", "[airy][xsf_tests]") {
    // Regression test: an overflow in acai (nn == -1) was returned as an underflow, i.e. a silent 0.
    // Like ZAIRY, airy now reports an overflow, so the wrapper returns NaN. |Ai(z)| ~ 1.6e302 and
    // |Ai'(z)| ~ 1.8e303 are both below DBL_MAX, so this overflow is early (as in AMOS).

    /*
    If you want to check the values with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    z = mpmath.mpc(-0.0, 129.94937845989423)
    print(mpmath.airyai(z), mpmath.airyai(z, derivative=1))
    */
    const std::complex<double> z(-0.0, 129.94937845989423);
    std::complex<double> ai, aip, bi, bip;
    xsf::airy(z, ai, aip, bi, bip);

    CAPTURE(z, ai, aip);
    REQUIRE(std::isnan(ai.real()));
    REQUIRE(std::isnan(aip.real()));
}
