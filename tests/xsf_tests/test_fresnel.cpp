#include "../testing_utils.h"

#include <xsf/cpu/fresnel.h>

#include <vector>

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

TEST_CASE("fcszo many zeros", "[fcszo][xsf_tests]") {
    // The deflation product prod(z - zo[i]) overflowed, so zeros 332 and up were NaN.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 40

    for kf, f in [(1, mpmath.fresnelc), (2, mpmath.fresnels)]:
        for n in [332, 1000]:
            # Newton from the asymptotic start value of the n-th zero
            psq = mpmath.sqrt(4 * n - 1) if kf == 1 else 2 * mpmath.sqrt(n)
            z0 = mpmath.mpc(
                psq - mpmath.log(mpmath.pi * psq) / (mpmath.pi**2 * psq**3),
                mpmath.log(mpmath.pi * psq) / (mpmath.pi * psq),
            )
            z = mpmath.findroot(f, z0)
            print(f"{kf}, {n}: {{{float(z.real)!r}, {float(z.imag)!r}}}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, int, std::complex<double>>;
    auto [kf, n, desired] = GENERATE(
        test_case{1, 332, {36.428026905677456, 0.04141828452596473}},
        test_case{1, 1000, {63.23765090489771, 0.026635556296725618}},
        test_case{2, 332, {36.44175000729524, 0.04140597780177087}},
        test_case{2, 1000, {63.2455570919721, 0.026632855870228636}}
    );

    std::vector<std::complex<double>> zo(n);
    xsf::fcszo(kf, n, zo.data());
    const double error = xsf::extended_relative_error(zo[n - 1], desired);
    const double rtol = 5e-15;
    CAPTURE(kf, n, zo[n - 1], desired, error, rtol);
    REQUIRE(error <= rtol);
}
