#include "../testing_utils.h"

#include <xsf/cpu/bessel.h>
#include <xsf/mdspan.h>

TEST_CASE("rctj large x", "[rctj][real][xsf_tests]") {
    // Regression test: the backward recurrence overflowed for x >~ 4476, so all orders were NaN.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 30

    for n, x in [(2, 5000), (5, 10000), (30, 20000), (5500, 5000)]:
        X = mpmath.mpf(x)
        # x * j_n(x) and its derivative
        j = lambda t: t * mpmath.sqrt(mpmath.pi / (2 * t)) * mpmath.besselj(n + mpmath.mpf(1) / 2, t)
        print(f"n={n}, x={x}: {float(j(X))!r} {float(mpmath.diff(j, X))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, double, double, double>;
    auto [n, x, desired_rj, desired_dj] = GENERATE(
        test_case{2, 5000.0, 0.9878735191670958, -0.15526114887616732},
        test_case{5, 10000.0, 0.9516959470409938, -0.30704225460358353},
        test_case{30, 20000.0, -0.600733011725298, -0.7994495137989018},
        test_case{5500, 5000.0, 2.2365323921212067e-65, 1.026728990451556e-65}
    );

    std::vector<double> rj(n + 1), dj(n + 1);
    xsf::rctj(x, xsf::cxx::mdspan(rj.data(), rj.size()), xsf::cxx::mdspan(dj.data(), dj.size()));

    const double error_rj = xsf::extended_relative_error(rj[n], desired_rj);
    const double error_dj = xsf::extended_relative_error(dj[n], desired_dj);
    const double rtol = 1e-12;
    CAPTURE(n, x, rj[n], desired_rj, error_rj, dj[n], desired_dj, error_dj, rtol);
    REQUIRE(error_rj <= rtol);
    REQUIRE(error_dj <= rtol);
}

TEST_CASE("rctj high orders", "[rctj][real][xsf_tests]") {
    // Regression test: orders above msta1(x, 200) were left unset, and the ones below it were inaccurate.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 30

    for k, x in [(105, 1), (45, "1e-3"), (176, 10)]:
        X = mpmath.mpf(x)
        # x * j_k(x) and its derivative
        j = lambda t: t * mpmath.sqrt(mpmath.pi / (2 * t)) * mpmath.besselj(k + mpmath.mpf(1) / 2, t)
        print(f"k={k}, x={x}: {float(j(X))!r} {float(mpmath.diff(j, X))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, int, double, double, double>;
    auto [n, k, x, desired_rj, desired_dj] = GENERATE(
        test_case{300, 105, 1.0, 1.9599728528102047e-201, 2.0774792044645683e-199},
        test_case{120, 45, 1e-3, 3.1130384250545827e-209, 1.4319976751903727e-204},
        test_case{300, 176, 10.0, 3.0545044405010714e-198, 5.397861821994432e-197},
        test_case{300, 300, 1.0, 0.0, 0.0} // underflows
    );

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> rj(n + 1, nan), dj(n + 1, nan);
    xsf::rctj(x, xsf::cxx::mdspan(rj.data(), rj.size()), xsf::cxx::mdspan(dj.data(), dj.size()));

    const double error_rj = xsf::extended_relative_error(rj[k], desired_rj);
    const double error_dj = xsf::extended_relative_error(dj[k], desired_dj);
    const double rtol = 1e-12;
    CAPTURE(n, k, x, rj[k], desired_rj, error_rj, dj[k], desired_dj, error_dj, rtol);
    REQUIRE(error_rj <= rtol);
    REQUIRE(error_dj <= rtol);
}
