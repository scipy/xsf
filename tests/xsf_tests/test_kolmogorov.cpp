#include "../testing_utils.h"

#include <xsf/cephes/kolmogorov.h>

TEST_CASE("cephes::kolmogorov and kolmogp just above the cutover", "[cephes::kolmogorov][xsf_tests]") {
    // Regression test for gh-451: for x > 0.82 the series dropped terms, so
    // kolmogp had a relative error up to 3e-8 and kolmogorov up to 1e-14.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50


    def sf(x):
        # DLMF 20.2.4: 1 - theta_4(0, exp(-2 x**2))
        return 1 - mpmath.jtheta(4, 0, mpmath.exp(-2 * x**2))


    for x in [0.8205, 0.83, 1.0]:
        x = mpmath.mpf(x)
        print(f"{float(x)!r}: {float(sf(x))!r}, {float(mpmath.diff(sf, x))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, desired_sf, desired_p] = GENERATE(
        test_case{0.8205, 0.5111775653643441, -1.5877555646537649},
        test_case{0.83, 0.4961909953505058, -1.5670710147480162},
        test_case{1.0, 0.2699996716773545, -1.0719485583569417}
    );

    const double out_sf = xsf::cephes::kolmogorov(x);
    const double out_p = xsf::cephes::kolmogp(x);
    const double error_sf = xsf::extended_relative_error(out_sf, desired_sf);
    const double error_p = xsf::extended_relative_error(out_p, desired_p);
    const double rtol = 5e-16;
    CAPTURE(x, out_sf, desired_sf, error_sf, out_p, desired_p, error_p, rtol);
    REQUIRE(error_sf <= rtol);
    REQUIRE(error_p <= rtol);
}
