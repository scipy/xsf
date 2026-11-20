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
