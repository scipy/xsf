#include "../testing_utils.h"

#include <cmath>
#include <complex>
#include <xsf/expint.h>

TEST_CASE("exp1, expi side of the cut for abs(z) >= 40", "[exp1][expi][xsf_tests]") {
    // The continued fraction branch ignored the sign of Im z on the cut, e.g. exp1(-40-0j)
    // had imaginary part -pi instead of +pi, and expi(40+0j) had 2 pi instead of 0.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    print(repr(float(mpmath.e1(-40).real)), repr(float(mpmath.ei(40).real)))
    */
    SET_FP_FORMAT()
    const double e1_re = -6039718263611242.0;
    const double ei_re = 6039718263611242.0;
    const double rtol = 5e-15;

    const double imag = GENERATE(0.0, -0.0);
    const double sign = std::signbit(imag) ? -1.0 : 1.0;

    // E1(-40 + 0j) = -Ei(40) - pi i, and E1(conj z) = conj E1(z).
    const std::complex<double> e1 = xsf::exp1(std::complex<double>(-40.0, imag));
    REQUIRE(xsf::extended_relative_error(e1.real(), e1_re) <= rtol);
    REQUIRE(std::abs(e1.imag() + sign * M_PI) <= 4 * 2.2e-16 * M_PI);

    // Ei(40 + 0j) is real.
    const std::complex<double> ei = xsf::expi(std::complex<double>(40.0, imag));
    REQUIRE(xsf::extended_relative_error(ei.real(), ei_re) <= rtol);
    REQUIRE(std::abs(ei.imag()) <= 4 * 2.2e-16 * M_PI);
}

TEST_CASE("exp1, expi overflow only above the true overflow threshold", "[exp1][expi][xsf_tests]") {
    // exp(x) overflows for x > 709.78 but Ei(x) only for x > 716.35, so these used to be inf/nan.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    print(repr(float(mpmath.ei(712).real)))
    e1 = mpmath.e1(mpmath.mpc(-712, 1))
    print(repr(float(e1.real)), repr(float(e1.imag)))
    e1 = mpmath.e1(mpmath.mpc(-716.65, 0.7853981633974483))
    print(repr(float(e1.real)), repr(float(e1.imag)))
    ei = mpmath.ei(mpmath.mpc(712, 1))
    print(repr(float(ei.real)), repr(float(ei.imag)))
    */
    SET_FP_FORMAT()
    const double rtol = 1e-14;

    const double ei = xsf::expi(712.0);
    REQUIRE(xsf::extended_relative_error(ei, 2.3216800841052115e306) <= rtol);

    const std::complex<double> e1 = xsf::exp1(std::complex<double>(-712.0, 1.0));
    REQUIRE(xsf::extended_relative_error(e1.real(), -1.2571543434669668e306) <= rtol);
    REQUIRE(xsf::extended_relative_error(e1.imag(), 1.9518582701792164e306) <= rtol);

    // Close to the threshold, where both components are near DBL_MAX.
    const std::complex<double> e1b = xsf::exp1(std::complex<double>(-716.65, 0.7853981633974483));
    REQUIRE(xsf::extended_relative_error(e1b.real(), -1.7076602595684751e308) <= rtol);
    REQUIRE(xsf::extended_relative_error(e1b.imag(), 1.7039161718563943e308) <= rtol);

    const std::complex<double> eic = xsf::expi(std::complex<double>(712.0, 1.0));
    REQUIRE(xsf::extended_relative_error(eic.real(), 1.2571543434669668e306) <= rtol);
    REQUIRE(xsf::extended_relative_error(eic.imag(), 1.9518582701792164e306) <= rtol);
}

TEST_CASE("exp1 near the negative real axis for abs(z) >= 40", "[exp1][xsf_tests]") {
    // Regression test: the continued fraction missed (part of) the 2 pi i jump across the cut.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 400

    for x, y in [(-40, -1e-10), (-45, 1e-17), (-100, -1e-50), (-710, 5e-310)]:
        e1 = mpmath.e1(mpmath.mpc(x, y))
        print(repr(float(e1.real)), repr(float(e1.imag)))
    */
    SET_FP_FORMAT()
    auto [x, y, re, im] = GENERATE(
        std::make_tuple(-40.0, -1e-10, -6039718263611242.0, -588460.0254998964),
        std::make_tuple(-45.0, 1e-17, -7.943916035704454e+17, 4.621578692518006),
        std::make_tuple(-100.0, -1e-50, -2.71555274485388e+41, 3.1415926509016763),
        std::make_tuple(-710.0, 5e-310, -3.1509156882062014e+305, -3.1414353300147115)
    );
    // exp(710) can be off by about 710 eps on other platforms
    const double rtol = 2e-13;

    const std::complex<double> e1 = xsf::exp1(std::complex<double>(x, y));
    CAPTURE(x, y, e1);
    REQUIRE(xsf::extended_relative_error(e1.real(), re) <= rtol);
    REQUIRE(xsf::extended_relative_error(e1.imag(), im) <= rtol);
}
