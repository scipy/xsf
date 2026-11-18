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

TEST_CASE("specfun::cyzo many real zeros", "[cyzo][xsf_tests]") {
    // Regression test: the deflation product prod(z - zo[i]) overflowed from about the 140th
    // zero on, so the 140th zero was off by ~4e-6 and later zeros were NaN.

    /*
    import mpmath as mp

    mp.mp.dps = 40
    for kf in range(3):
        if kf == 0:
            x140, x150 = mp.besselyzero(0, 140), mp.besselyzero(0, 150)
            v = mp.bessely(1, x150)
        elif kf == 1:
            x140, x150 = mp.besselyzero(1, 140), mp.besselyzero(1, 150)
            v = mp.bessely(0, x150)
        else:
            x140, x150 = mp.besselyzero(1, 140, 1), mp.besselyzero(1, 150, 1)
            v = mp.bessely(1, x150)
        print(kf, repr(float(x140)), repr(float(x150)), repr(float(v)))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, double, double, double>;
    auto [kf, desired_zero140, desired_zero150, desired_value150] = GENERATE(
        test_case{0, 437.4670627473985, 468.882970138688, 0.0368474938665893},
        test_case{1, 439.03671919843214, 470.4527027720758, -0.03678592741791805},
        test_case{2, 440.6063837622912, 472.022442473055, -0.03672475095937137}
    );

    std::complex<double> zo[150] = {};
    std::complex<double> zv[150] = {};
    xsf::specfun::cyzo(150, kf, 1, zo, zv);

    const double error140 = xsf::extended_relative_error(zo[139].real(), desired_zero140);
    const double error150 = xsf::extended_relative_error(zo[149].real(), desired_zero150);
    const double error_value = xsf::extended_relative_error(zv[149].real(), desired_value150);
    CAPTURE(kf, zo[139], zo[149], zv[149], error140, error150, error_value);
    REQUIRE(error140 <= 1e-14);
    REQUIRE(error150 <= 1e-14);
    REQUIRE(error_value <= 1e-13);
}
