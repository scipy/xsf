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

TEST_CASE("sici complex Ci for large Re(z) and at Im(z) = -0", "[sici][xsf_tests]") {
    // Regression test: Ci lost its tiny imaginary part for huge Re(z), and
    // had Im(Ci) = +pi instead of -pi on the negative real axis at Im(z) = -0.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 80

    z = complex(2.981344618416097e61, 55.61111989001531)
    print(f"Ci({z!r}) = {complex(mpmath.ci(z))!r}")
    # Ci(conj(z)) = conj(Ci(z)), also for Im(z) = -0
    print(f"Ci((-1-0j)) = {complex(mpmath.ci(-1 + 0j)).conjugate()!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [x, y, desired_re, desired_im] = GENERATE(
        test_case{2.981344618416097e61, 55.61111989001531, -7.448138578702911e-39, -2.2580510795720128e-38},
        test_case{-1.0, -0.0, 0.33740392290096816, -3.141592653589793},
        test_case{-std::numeric_limits<double>::infinity(), -0.0, 0.0, -3.141592653589793}
    );

    std::complex<double> si, ci;
    xsf::sici(std::complex<double>(x, y), si, ci);
    const double error_re = xsf::extended_relative_error(ci.real(), desired_re);
    const double error_im = xsf::extended_relative_error(ci.imag(), desired_im);
    const double rtol = 1e-14;
    CAPTURE(x, y, ci, desired_re, desired_im, error_re, error_im, rtol);
    REQUIRE(error_re <= rtol);
    REQUIRE(error_im <= rtol);
}
