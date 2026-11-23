#define SP_SPECFUN_ERROR
#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/cpu/bessel.h>

namespace {
sf_error_t last_error = SF_ERROR_OK;
}

void xsf::set_error(const char *, sf_error_t code, const char *, ...) { last_error = code; }

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

TEST_CASE("cyl_bessel_k and cyl_bessel_ke huge order", "[cyl_bessel_k][cyl_bessel_ke][xsf_tests]") {
    // Regression test for gh-499: these returned NaN instead of inf.
    SET_FP_FORMAT()
    const double inf = std::numeric_limits<double>::infinity();
    auto [v, x] = GENERATE(std::tuple{0x1p30, 1.0}, std::tuple{-1e12, 1.0}, std::tuple{1e300, 1e200});
    CAPTURE(v, x);
    REQUIRE(xsf::cyl_bessel_k(v, x) == inf);
    REQUIRE(xsf::cyl_bessel_ke(v, x) == inf);
    REQUIRE(xsf::cyl_bessel_k(v, std::complex<double>(x)) == std::complex<double>(inf, 0.0));
}

TEST_CASE("kve huge order with x > v", "[cyl_bessel_k][xsf_tests]") {
    // Regression test for gh-499: kve(v, x) overflows for x > v too (log kve = 3.5e8 here).
    SET_FP_FORMAT()
    REQUIRE(xsf::cyl_bessel_ke(0x1p30, 0x1.8p30) == std::numeric_limits<double>::infinity());
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

    last_error = SF_ERROR_OK;
    const double k = xsf::cyl_bessel_k(v, x);
    const double ke = xsf::cyl_bessel_ke(v, x);
    REQUIRE(last_error == SF_ERROR_OK);
    const double k_error = xsf::extended_relative_error(k, desired_k);
    const double ke_error = xsf::extended_relative_error(ke, desired_ke);
    // AMOS itself is only accurate to a few 1e-14 here, and other libms differ
    const double rtol = 5e-14;
    CAPTURE(v, x, k, ke, desired_k, desired_ke, k_error, ke_error, rtol);
    REQUIRE(k_error <= rtol);
    REQUIRE(ke_error <= rtol);
}

TEST_CASE("cyl_bessel_k and cyl_bessel_ke overflow error", "[cyl_bessel_k][cyl_bessel_ke][xsf_tests]") {
    // Regression test for gh-372: the overflow error is still reported where the result is inf.
    last_error = SF_ERROR_OK;
    REQUIRE(xsf::cyl_bessel_k(1.5, 1e-206) == std::numeric_limits<double>::infinity());
    REQUIRE(last_error == SF_ERROR_OVERFLOW);
    last_error = SF_ERROR_OK;
    REQUIRE(xsf::cyl_bessel_ke(1.5, 1e-206) == std::numeric_limits<double>::infinity());
    REQUIRE(last_error == SF_ERROR_OVERFLOW);
}

TEST_CASE("cyl_bessel_k near the underflow threshold", "[cyl_bessel_k][xsf_tests]") {
    // Regression test: these returned 0 although K is a normal double.
    /*
    Reference values from mpmath:

    from mpmath import mp, besselk, mpc
    mp.dps = 30
    for v, z in [(0, 699), (50, 700), (0.5, mpc(702, -50)), (0, mpc(700, 1e5))]:
        print(repr(complex(besselk(v, z))))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, std::complex<double>, std::complex<double>>;
    auto [v, z, desired] = GENERATE(
        test_case{0.0, {699.0, 0.0}, {1.2702841880327418e-305, 0.0}},
        test_case{50.0, {700.0, 0.0}, {2.7793358770120587e-305, 0.0}},
        test_case{0.5, {702.0, -50.0}, {6.138086236570686e-307, -1.4367386773346883e-307}},
        test_case{0.0, {700.0, 1e5}, {-2.86945346472656e-307, 2.6525609459688144e-307}}
    );

    last_error = SF_ERROR_OK;
    const std::complex<double> w = xsf::cyl_bessel_k(v, z);
    const double error = xsf::extended_relative_error(w, desired);
    const double rtol = 5e-14;
    CAPTURE(v, z, w, desired, error, rtol);
    REQUIRE(error <= rtol);
    REQUIRE(last_error != SF_ERROR_UNDERFLOW);
    if (z.imag() == 0) {
        const double k = xsf::cyl_bessel_k(v, z.real());
        const double k_error = xsf::extended_relative_error(k, desired.real());
        CAPTURE(k, k_error);
        REQUIRE(k_error <= rtol);
    }
}
