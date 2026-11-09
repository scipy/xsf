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
