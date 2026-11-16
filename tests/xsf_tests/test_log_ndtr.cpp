#include "../testing_utils.h"

#include <xsf/stats.h>

TEST_CASE("log_ndtr complex z with small ndtr(-z)", "[log_ndtr][xsf_tests]") {
    // Regression test: for Re(z) > 0 and small ndtr(-z), log(ndtr(z)) = log(1 - ndtr(-z))
    // was computed as a sum that cancels, with relative errors up to 1e-6.

    /*
    import mpmath

    mpmath.mp.dps = 50

    for z in [5.901265878562874 + 0.4169925459625219j, 4.5 - 1.25j, 6.25 - 0.75j]:
        # log(ndtr(z)) = log1p(-ndtr(-z)), with ndtr(-z) = erfc(z / sqrt(2)) / 2
        w = mpmath.log1p(-mpmath.erfc(z / mpmath.sqrt(2)) / 2)
        print(f"{z!r}: {float(w.real)!r}, {float(w.imag)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<std::complex<double>, std::complex<double>>;
    auto [z, desired] = GENERATE(
        test_case{{5.901265878562874, 0.4169925459625219}, {1.6048642384381675e-09, 1.130791271554178e-09}},
        test_case{{4.5, -1.25}, {-6.615497872514756e-06, 2.8530732616772225e-06}},
        test_case{{6.25, -0.75}, {-2.4067183330437477e-11, 2.691130913270477e-10}}
    );

    const std::complex<double> out = xsf::log_ndtr(z);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-14;
    CAPTURE(z, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
