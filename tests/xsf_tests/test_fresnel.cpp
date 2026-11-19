#include "../testing_utils.h"

#include <xsf/cpu/fresnel.h>

TEST_CASE("fresnel huge x", "[fresnel][xsf_tests]") {
    // Regression test for gh-386: this returned NaN for |x| >= 2^512, where x * x overflows.
    // mpmath.fresnels(x) and mpmath.fresnelc(x) round to +-0.5 here.
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(test_case{0x1p512, 0.5}, test_case{1e155, 0.5}, test_case{-1e200, -0.5});

    double s, c;
    xsf::fresnel(x, s, c);
    CAPTURE(x, s, c, desired);
    REQUIRE(s == desired);
    REQUIRE(c == desired);
}

TEST_CASE("modified_fresnel_plus/minus K large x", "[modified_fresnel][xsf_tests]") {
    // Regression test for gh-453: K+-(x) lost digits like eps * x^2 for x >= 5.5, and was NaN once x * x
    // overflowed.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for x in [1e4, 1e9, 1e200]:
        # K+(x) = exp(-i x^2) erfc(exp(-i pi/4) x) / 2, DLMF 7.5.10; K-(x) = conj(K+(x))
        with mpmath.workdps(450):
            k = mpmath.exp(-1j * mpmath.mpf(x) ** 2) * mpmath.erfc(mpmath.expjpi(-0.25) * x) / 2
        print(f"K+({x!r}) = complex({float(k.real)!r}, {float(k.imag)!r})")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, std::complex<double>>;
    auto [x, desired] = GENERATE(
        test_case{1e4, {1.9947114119807202e-05, 1.9947113920336063e-05}},
        test_case{1e9, {1.9947114020071633e-10, 1.9947114020071633e-10}},
        test_case{1e200, {1.9947114020071636e-201, 1.9947114020071636e-201}}
    );

    std::complex<double> fp, kp, fm, km;
    xsf::modified_fresnel_plus(x, fp, kp);
    xsf::modified_fresnel_minus(x, fm, km);
    const double error_p = xsf::extended_relative_error(kp, desired);
    const double error_m = xsf::extended_relative_error(km, std::conj(desired));
    const double rtol = 1e-15;
    CAPTURE(x, kp, km, desired, error_p, error_m, rtol);
    REQUIRE(error_p <= rtol);
    REQUIRE(error_m <= rtol);
}
