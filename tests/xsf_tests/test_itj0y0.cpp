#include "../testing_utils.h"
#include <tuple>
#include <xsf/bessel.h>

TEST_CASE("itj0y0", "[itj0y0][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double>;
    // Reference values were computed with the Python library mpmath.
    auto [x, ref_j0, ref_y0, rtol] = GENERATE(
        test_case{10, 1.0670113039567368, 0.24129031832266684, 1e-11},
        test_case{20, 1.0583788214211278, -0.1682159767721503, 1e-7},
        test_case{30, 0.8842490888254749, 0.08822971194803665, 1e-11},
        test_case{40, 1.1257761503599915, -0.008932466273627416, 1e-11},
        test_case{50, 0.9014121225818346, -0.05481407100006349, 1e-11}
    );
    double j0int, y0int;
    xsf::it1j0y0(x, j0int, y0int);

    const auto rel_error_j0 = xsf::extended_relative_error(j0int, ref_j0);
    const auto rel_error_y0 = xsf::extended_relative_error(y0int, ref_y0);

    CAPTURE(x, j0int, ref_j0, rtol, rel_error_j0);
    CAPTURE(x, y0int, ref_y0, rtol, rel_error_y0);
    REQUIRE(rel_error_j0 <= rtol);
    REQUIRE(rel_error_y0 <= rtol);
}
