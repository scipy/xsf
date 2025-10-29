#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/bessel.h>

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
