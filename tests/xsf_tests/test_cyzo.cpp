#include "../testing_utils.h"

#include <xsf/cpu/specfun/specfun.h>

TEST_CASE("specfun::cyzo real zeros with |z| > 12", "[cyzo][xsf_tests]") {
    // Regression test: the asymptotic coefficient tables in cy01 (used for |z| > 12) lacked
    // their exponents, so zeros of Y0, Y1 and Y1' beyond |z| = 12 were off by up to ~1e-9.

    /*
    import mpmath as mp

    mp.mp.dps = 40
    for kf in range(3):
        if kf == 0:
            x = mp.besselyzero(0, 5)
            v = mp.bessely(1, x)
        elif kf == 1:
            x = mp.besselyzero(1, 5)
            v = mp.bessely(0, x)
        else:
            x = mp.besselyzero(1, 5, derivative=1)
            v = mp.bessely(1, x)
        print(kf, repr(float(x)), repr(float(v)))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, double, double>;
    auto [kf, desired_zero, desired_value] = GENERATE(
        test_case{0, 13.361097473872764, -0.21835829659767134}, test_case{1, 14.897442128336726, 0.2065471103565926},
        test_case{2, 16.44005800729328, 0.19682760568385493}
    );

    std::complex<double> zo[5] = {};
    std::complex<double> zv[5] = {};
    xsf::specfun::cyzo(5, kf, 1, zo, zv);

    const double error = xsf::extended_relative_error(zo[4].real(), desired_zero);
    const double error_value = xsf::extended_relative_error(zv[4].real(), desired_value);
    const double rtol = 1e-12;
    CAPTURE(kf, zo[4], desired_zero, error, zv[4], desired_value, error_value, rtol);
    REQUIRE(error <= rtol);
    REQUIRE(error_value <= 5e-12);
}
