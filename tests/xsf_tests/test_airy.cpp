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

TEST_CASE("airy Ai Taylor", "[airy][xsf_tests]") {
    // Regression test: the Maclaurin series of Ai and Ai' cancel for 0.75 <= x < 2.09, with errors of up to ~100 eps.

    /*
    If you want to check the values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for x in (2.0799248120300753, 2.0631328320802003):
        x = mpmath.mpf(x)
        print(float(mpmath.airyai(x)), float(mpmath.airyai(x, derivative=1)))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, ai_desired, aip_desired] = GENERATE(
        test_case{2.0799248120300753, 0.03089798121053051, -0.04773240803820107},
        test_case{2.0631328320802003, 0.03170861742551274, -0.04882122346344751}
    );

    double ai, aip, bi, bip;
    xsf::airy(x, ai, aip, bi, bip);
    const double rtol = 2e-15;
    const double ai_error = xsf::extended_relative_error(ai, ai_desired);
    const double aip_error = xsf::extended_relative_error(aip, aip_desired);
    CAPTURE(x, ai, ai_desired, ai_error, aip, aip_desired, aip_error, rtol);
    REQUIRE(ai_error <= rtol);
    REQUIRE(aip_error <= rtol);
}
