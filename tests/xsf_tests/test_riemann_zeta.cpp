#include "../testing_utils.h"

#include <xsf/zeta.h>

TEST_CASE("riemann_zeta complex z near 0", "[riemann_zeta][xsf_tests]") {
    // Regression test: the reflection formula evaluated zeta(1 - z) next to its pole,
    // losing accuracy for non-real abs(z) < 1/2.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    from mpmath import mp, mpc, zeta

    mp.dps = 50
    for z in [1e-6 - 2e-6j, 1e-5 + 1e-5j, -1e-3 + 3e-3j]:
        w = zeta(mpc(z))
        print(f"{z!r}: ({float(w.real)!r}, {float(w.imag)!r})")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<std::complex<double>, std::complex<double>>;
    auto [z, desired] = GENERATE(
        test_case{{1e-6, -2e-6}, {-0.5000009189355237, 1.8378810791202557e-06}},
        test_case{{1e-5, 1e-5}, {-0.50000918938533, -9.189585969693889e-06}},
        test_case{{-1e-3, 3e-3}, {-0.4990730620890677, -0.002750778612088271}}
    );

    const auto out = xsf::riemann_zeta(z);
    // The imaginary part is checked separately: it is much smaller than the real part.
    const double re_error = xsf::extended_relative_error(out.real(), desired.real());
    const double im_error = xsf::extended_relative_error(out.imag(), desired.imag());
    const double rtol = 5e-15;
    CAPTURE(z, out, desired, re_error, im_error, rtol);
    REQUIRE(re_error <= rtol);
    REQUIRE(im_error <= rtol);
}

TEST_CASE("riemann_zeta(z) for 30 < |Im z| <= 50 and near the zeros of 1 - 2^(1 - z)", "[riemann_zeta][xsf_tests]") {
    // Borwein's 50-term sum is inaccurate for |Im z| > 30, and dividing by 1 - 2^(1 - z) loses
    // digits near its zeros 1 + 2 pi i k / ln(2). Also near the pole, where Re zeta is finite
    // while Im zeta is about -1e307, and Im zeta is much smaller than Re zeta at large Re z.
    // Reference values, with mpmath.mp.dps = 50:
    //     w = mpmath.zeta(mpmath.mpc(re, im)); print(float(w.real), float(w.imag))
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double, double>;
    auto [re, im, desired_re, desired_im, rtol] = GENERATE(
        test_case{0.5, 48.0, -0.006058168095158761, -0.00536975820007769, 1e-11},
        test_case{-0.5, 47.0, -5.485325059905793, -7.224137203791469, 1e-12},
        test_case{33.0, 40.0, 0.9999999999006589, -6.069455759733579e-11, 1e-12},
        test_case{27.0, 31.726520992790356, 0.9999999925492941, 3.845617282366759e-14, 1e-9},
        test_case{1.0, 9.064720283654388, 1.3465795428363172, 0.1098831367962695, 1e-12},
        test_case{1.0, 1e-307, 0.5772156649015329, -1.0000000000000001e+307, 1e-12}
    );

    const std::complex<double> out = xsf::riemann_zeta(std::complex<double>(re, im));
    const double error_re = xsf::extended_relative_error(out.real(), desired_re);
    const double error_im = xsf::extended_relative_error(out.imag(), desired_im);
    CAPTURE(re, im, out, desired_re, desired_im, error_re, error_im, rtol);
    REQUIRE(error_re <= rtol);
    REQUIRE(error_im <= rtol);
}
