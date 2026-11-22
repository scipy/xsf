#include "../testing_utils.h"

#include <xsf/wright_bessel.h>

TEST_CASE("log_wright_bessel result close to 0", "[log_wright_bessel][xsf_tests]") {
    // Test that log_wright_bessel keeps the precision when Phi is close to 1, i.e.
    // when the result is close to 0. The last three cases have b one ulp below 1,
    // where log(Gamma(b)) is close to 0 as well.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 60

    def log_wright_bessel(a, b, x):
        a, b, x = map(mpmath.mpf, (a, b, x))
        t0 = mpmath.rgamma(b)
        # sum over k >= 1, so that log1p keeps the small terms
        tail = mpmath.nsum(lambda k: x**k / mpmath.factorial(k) * mpmath.rgamma(a * k + b), [1, mpmath.inf])
        return mpmath.log(t0) + mpmath.log1p(tail / t0)

    b1 = 0.9999999999999999
    for a, b, x in [
        (20, 1, 2.5), (100, 1, 1.5), (10, 2, 1e-3), (1e-4, 1, 1e-10),
        (10, b1, 3.63e-10), (100, b1, 2.5), (1e-4, b1, 1e-20),
    ]:
        print(f"log_wright_bessel({a}, {b}, {x}) = {float(log_wright_bessel(a, b, x))!r}")
    */
    SET_FP_FORMAT()
    static constexpr double b1 = 0.9999999999999999; // the largest double below 1
    using test_case = std::tuple<double, double, double, double>;
    auto [a, b, x, desired] = GENERATE(
        test_case{20, 1, 2.5, 1.0275794058280412e-18}, test_case{100, 1, 1.5, 1.6072654321882003e-158},
        test_case{10, 2, 1e-3, 2.5052108385127925e-11}, test_case{1e-4, 1, 1e-10, 1.0000577150076675e-10},
        test_case{10, b1, 3.63e-10, 3.594925664826874e-17}, test_case{100, b1, 2.5, -6.408381213480008e-17},
        test_case{1e-4, b1, 1e-20, -6.407381155765e-17}
    );

    const double out = xsf::log_wright_bessel(a, b, x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-12;
    CAPTURE(a, b, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
