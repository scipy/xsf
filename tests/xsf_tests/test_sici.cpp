#include "../testing_utils.h"

#include <cmath>
#include <complex>
#include <xsf/sici.h>

TEST_CASE("sici imaginary-axis branches", "[sici][xsf_tests]") {
    const double real = GENERATE(0.0, -0.0);
    const double imag = GENERATE(-10.0, -1.0, 1.0, 10.0);
    const std::complex<double> z(real, imag);
    std::complex<double> si, ci;
    xsf::sici(z, si, ci);

    // Si(iy) = i Shi(y), Ci(iy) = Chi(|y|) + sign(y) i pi/2.
    // Use the independent real-valued implementation as the reference.
    double shi, chi;
    xsf::shichi(std::abs(imag), shi, chi);
    CAPTURE(real, imag, si, ci);
    REQUIRE(std::abs(si.real()) <= 1e-13);
    REQUIRE(xsf::extended_relative_error(si.imag(), std::copysign(shi, imag)) <= 1e-13);
    REQUIRE(xsf::extended_relative_error(ci.real(), chi) <= 1e-13);
    REQUIRE(std::abs(ci.imag() - std::copysign(M_PI_2, imag)) <= 1e-13);
}
