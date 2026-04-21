#include "../testing_utils.h"
#include <xsf/cephes/j0.h>

TEST_CASE("j0 right tail gh-large-input", "[j0][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    // Reference values computed with mpmath with 1000 digts of precision
    auto [x, ref] = GENERATE(
        test_case{1e15, 6.156638646885021e-09}, test_case{1e12, 1.0167125050040682e-07},
        test_case{1e13, 1.1926484739665653e-07}, test_case{1e8, 3.206029534041208e-05},
        test_case{1e20, 6.698009040703424e-12}
    );
    const double w = xsf::cephes::j0(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rel_error);
    REQUIRE(rel_error <= 1e-15);
}
