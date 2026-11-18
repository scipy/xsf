#include "../testing_utils.h"
#include <tuple>
#include <xsf/stats.h>

TEST_CASE("owens_t large h small a", "[owens_t][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
    // Reference values were computed with the Python library mpmath.
    auto [h, a, ref] = GENERATE(
        test_case{5.0, 1e-10, 5.931152735254122e-17}, test_case{6.0, 1e-14, 2.423926559560457e-23},
        test_case{8.0, 1e-3, 2.0155423762820997e-18}
    );
    constexpr double rtol = 1e-13;
    const auto out = xsf::owens_t(h, a);
    const auto rel_error = xsf::extended_relative_error(out, ref);
    CAPTURE(h, a, out, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
