#include "../testing_utils.h"

#include <complex>
#include <xsf/sici.h>

TEST_CASE("shichi on the negative real axis", "[shichi][xsf_tests]") {
    // The side of the branch cut follows the sign of Im(z), also for -0:
    // Chi(-x +- 0i) = Chi(x) +- i pi, Shi(-x) = -Shi(x).

    /*
    import mpmath

    mpmath.mp.dps = 50

    for x in [1.0, 10.0]:
        print(f"Shi({x}) = {float(mpmath.shi(x))!r}, Chi({x}) = {float(mpmath.chi(x))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, shi_x, chi_x] = GENERATE(
        test_case{1.0, 1.0572508753757286, 0.8378669409802082}, test_case{10.0, 1246.1144901994232, 1246.1144860424545}
    );
    const double imag = GENERATE(0.0, -0.0);

    std::complex<double> shi, chi;
    xsf::shichi(std::complex<double>(-x, imag), shi, chi);
    const double rtol = 5e-15;
    CAPTURE(x, imag, shi, chi);
    REQUIRE(xsf::extended_relative_error(shi.real(), -shi_x) <= rtol);
    REQUIRE(xsf::extended_relative_error(chi.real(), chi_x) <= rtol);
    REQUIRE(shi.imag() == 0.0);
    REQUIRE(chi.imag() == std::copysign(M_PI, imag));
}
