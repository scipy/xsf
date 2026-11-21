#include "../testing_utils.h"

#include <cmath>
#include <tuple>
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

TEST_CASE("airyzo fifth zeros", "[airyzo][xsf_tests]") {
    // Regression test: airyzo refined the zeros with airyb, whose power series cancels for
    // -8 <= x < 0, so the 5th zeros had relative errors of up to 1.6e-11.

    /*
    If you want to check the values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    k = 5
    a = mpmath.airyaizero(k)
    ap = mpmath.airyaizero(k, derivative=1)
    b = mpmath.airybizero(k)
    bp = mpmath.airybizero(k, derivative=1)
    print(float(a), float(ap), float(mpmath.airyai(ap)), float(mpmath.airyai(a, derivative=1)))
    print(float(b), float(bp), float(mpmath.airybi(bp)), float(mpmath.airybi(b, derivative=1)))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, std::array<double, 4>>;
    auto [kf, desired] = GENERATE(
        test_case{1, {-7.944133587120853, -7.37217725504777, 0.34230124441162385, 0.9473357094415678}},
        test_case{2, {-7.376762079367763, -7.940178689168579, -0.33602624013366006, 0.9299836385680267}}
    );

    double xa[5], xb[5], xc[5], xd[5];
    xsf::airyzo(5, kf, xa, xb, xc, xd);
    const double out[4] = {xa[4], xb[4], xc[4], xd[4]};
    for (int j = 0; j < 4; ++j) {
        const double error = xsf::extended_relative_error(out[j], desired[j]);
        const double rtol = 1e-14;
        CAPTURE(kf, j, out[j], desired[j], error, rtol);
        REQUIRE(error <= rtol);
    }
}

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
