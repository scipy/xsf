#include "../testing_utils.h"

#include <xsf/cpu/bessel.h>

TEST_CASE("besselpoly large a", "[besselpoly][xsf_tests]") {
    // Regression test for gh-452: the power series in -a**2 lost digits for
    // a > 3 and diverged for larger a.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for a, lmb, nu in [(7.0, 0.0, 0.0), (30.0, 1.0, 0.5)]:
        f = lambda t: t**lmb * mpmath.besselj(nu, 2 * a * t)
        res = mpmath.quad(f, mpmath.linspace(0, 1, 50))
        print(f"besselpoly({a!r}, {lmb!r}, {nu!r}) = {float(res)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [a, lmb, nu, desired] =
        GENERATE(test_case{7.0, 0.0, 0.0, 0.08004841666588781}, test_case{30.0, 1.0, 0.5, 0.0017697205282620003});

    const double out = xsf::besselpoly(a, lmb, nu);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-14;
    CAPTURE(a, lmb, nu, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("besselpoly large integer nu", "[besselpoly][xsf_tests]") {
    // Regression test: cephes::jv overflowed in its backward recurrence for
    // integer orders where J is tiny, so the Neumann series returned NaN.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for a, lmb, nu in [(42.23, 0.0, 484.0), (43.75, 0.0, 490.0)]:
        c = mpmath.mpf(lmb + nu + 1) / 2
        res = mpmath.hyp1f2(c, nu + 1, c + 1, -mpmath.mpf(a)**2)
        res *= mpmath.mpf(a)**nu / (2 * c * mpmath.gamma(nu + 1))
        print(f"besselpoly({a!r}, {lmb!r}, {nu!r}) = {float(res)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [a, lmb, nu, desired] = GENERATE(
        test_case{42.23, 0.0, 484.0, 3.23801322587164e-309}, test_case{43.75, 0.0, 490.0, 3.6250468855274235e-308}
    );

    const double out = xsf::besselpoly(a, lmb, nu);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-12;
    CAPTURE(a, lmb, nu, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("besselpoly large nu Neumann series", "[besselpoly][xsf_tests]") {
    // Follow-up to gh-485: the Neumann series was normalized with cephes::jv, which
    // loses digits or underflows to 0 for large nu.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for a, lmb, nu in [(60.0, 0.0, 244.5), (40.0, 0.0, 400.5)]:
        c = mpmath.mpf(lmb + nu + 1) / 2
        res = mpmath.hyp1f2(c, nu + 1, c + 1, -mpmath.mpf(a)**2)
        res *= mpmath.mpf(a)**nu / (2 * c * mpmath.gamma(nu + 1))
        print(f"besselpoly({a!r}, {lmb!r}, {nu!r}) = {float(res)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [a, lmb, nu, desired] = GENERATE(
        test_case{60.0, 0.0, 244.5, 3.230211041251911e-54}, test_case{40.0, 0.0, 400.5, 1.5235998731386507e-233}
    );

    const double out = xsf::besselpoly(a, lmb, nu);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-13;
    CAPTURE(a, lmb, nu, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
