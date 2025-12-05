#include "../testing_utils.h"
#include <tuple>
#include <xsf/legendre.h>

TEST_CASE("assoc_legendre_p scipy/gh-23101", "[assoc_legendre_p][xsf_tests]") {
    using test_case = std::tuple<int, int, double, double, double>;
    // Reference values were computed with the Python library mpmath.
    auto [n, m, z, ref, rtol] = GENERATE(
        test_case{1, 0, -1.0, -1.0, 1e-15}, test_case{1, 0, std::nextafter(-1.0, -2.0), -1.0000000000000002, 1e-15},
        test_case{1, 0, std::nextafter(-1.0, 0.0), -0.9999999999999999, 1e-15}, test_case{2, 0, -1.0, 1.0, 1e-15},
        test_case{2, 0, std::nextafter(-1.0, -2.0), 1.0000000000000007, 1e-15},
        test_case{2, 0, std::nextafter(-1.0, 0.0), 0.9999999999999997, 1e-15}, test_case{3, 0, -1.0, -1.0, 1e-15},
        test_case{3, 0, std::nextafter(-1.0, -2.0), -1.0000000000000013, 1e-15},
        test_case{3, 0, std::nextafter(-1.0, 0.0), -0.9999999999999993, 1e-15}, test_case{4, 0, -1.0, 1.0, 1e-15},
        test_case{4, 0, std::nextafter(-1.0, -2.0), 1.0000000000000022, 4e-15},
        test_case{4, 0, std::nextafter(-1.0, 0.0), 0.9999999999999989, 1e-15}
    );
    const double w = xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, 2);
    const auto rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(n, m, z, w, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
