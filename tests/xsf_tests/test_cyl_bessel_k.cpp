#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/cpu/bessel.h>

TEST_CASE("cyl_bessel_k gh-46", "[cyl_bessel_k][xsf_tests]") {
    using test_case = std::tuple<double, std::complex<double>, std::complex<double>, double>;
    using std::complex;
    // Reference values were computed with the Python library mpmath.
    auto [v, z, ref, rtol] = GENERATE(
        test_case{0.0, complex{680.0, -1000.0}, complex{1.901684871999608e-298, 1.713412341479591e-297}, 1e-13},
        test_case{0.0, complex{680.0, -680.0}, complex{-4.553730032944803e-298, 1.878727010109855e-297}, 1e-13},
        test_case{0.0, complex{25.0, 100.0}, complex{1.699267100365868e-12, -2.234890030902166e-13}, 5e-16}
    );
    const complex w = xsf::cyl_bessel_k(v, z);
    const auto rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(v, z, w, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("cyl_bessel_k and cyl_bessel_ke near DBL_MAX", "[cyl_bessel_k][cyl_bessel_ke][xsf_tests]") {
    // Regression test for gh-372: these returned inf.
    /*
    Reference values from mpmath:

    from mpmath import mp, besselk, exp
    mp.dps = 30
    for v, x in [(1.5, 1e-205), (-70.0, 0.002)]:
        print(repr(float(besselk(v, x))), repr(float(exp(x) * besselk(v, x))))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [v, x, desired_k, desired_ke] = GENERATE(
        test_case{1.5, 1e-205, 3.963327297606011e+307, 3.963327297606011e+307},
        test_case{-70.0, 0.002, 8.556122497405277e+307, 8.573251866058952e+307}
    );

    const double k = xsf::cyl_bessel_k(v, x);
    const double ke = xsf::cyl_bessel_ke(v, x);
    const double k_error = xsf::extended_relative_error(k, desired_k);
    const double ke_error = xsf::extended_relative_error(ke, desired_ke);
    const double rtol = 5e-15;
    CAPTURE(v, x, k, ke, desired_k, desired_ke, k_error, ke_error, rtol);
    REQUIRE(k_error <= rtol);
    REQUIRE(ke_error <= rtol);
}
