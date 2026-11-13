#include "../testing_utils.h"
#include <cmath>
#include <tuple>
#include <xsf/cpu/airy.h>

TEST_CASE("airy real x < -10 phase", "[airy][xsf_tests]") {
    // Regression test for gh-532 (https://github.com/scipy/xsf/issues/532):
    // the phase of the oscillation lost digits, so the error relative to the
    // amplitude was ~1e-12 at x = -1e3, ~1e-9 at x = -1e5 and NaN at x = -1e7.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 60

    for x in [-1e3, -1e5, -1e7]:
        x = mpmath.mpf(x)
        print(
            f"{float(x)!r}: {float(mpmath.airyai(x))!r}, "
            f"{float(mpmath.airyai(x, derivative=1))!r}, {float(mpmath.airybi(x))!r}, "
            f"{float(mpmath.airybi(x, derivative=1))!r}"
        )
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double, double>;
    auto [x, ai, aip, bi, bip] = GENERATE(
        test_case{-1000.0, 0.05597189577301992, 2.6330710195241287, -0.08326457411708063, 1.769965940135989},
        test_case{-100000.0, -0.013152978737498166, -9.130078851932488, 0.02887184428505844, -4.1593370104464045},
        test_case{-10000000.0, 0.005418514944210602, -26.701707811976565, 0.008443821410258536, 17.134848759576794}
    );

    double out_ai, out_aip, out_bi, out_bip;
    xsf::airy(x, out_ai, out_aip, out_bi, out_bip);

    // The error is measured relative to the amplitude, as the outputs have zeros.
    const double amp = std::hypot(ai, bi);
    const double amp_p = std::hypot(aip, bip);
    const double rtol = 2e-15;
    const double error = (std::abs(out_ai - ai) + std::abs(out_bi - bi)) / amp +
                         (std::abs(out_aip - aip) + std::abs(out_bip - bip)) / amp_p;
    CAPTURE(x, out_ai, out_aip, out_bi, out_bip, error, rtol);
    REQUIRE(error <= rtol);
}
