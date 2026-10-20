#include "../testing_utils.h"
#include <tuple>
#include <xsf/cpu/struve.h>

TEST_CASE("itstruve0 large x", "[itstruve0][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
    // Reference values were computed with the Python library mpmath.
    auto [x, ref, rtol] = GENERATE(
        test_case{31.0, 2.946198073862365, 1e-13}, test_case{40.0, 3.148417627447186, 1e-13},
        test_case{100.0, 3.720914252854685, 1e-13}, test_case{300.0, 4.47323362427753, 1e-13}
    );
    const auto out = xsf::itstruve0(x);
    const auto rel_error = xsf::extended_relative_error(out, ref);
    CAPTURE(x, out, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
