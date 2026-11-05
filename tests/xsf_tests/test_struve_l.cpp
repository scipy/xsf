#include "../testing_utils.h"

#include <xsf/cpu/struve.h>

TEST_CASE("struve_l large-z expansion", "[struve_l][xsf_tests]") {
    // Regression test for gh-359: the large-z expansion added I_v instead of I_{-v}.
    // Reference values: float(mpmath.struvel(v, x)) with mpmath.mp.dps = 50.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [v, x, desired] =
        GENERATE(test_case{-15.5000001, 15, 162.7429503008563}, test_case{0.5, 12.5, 30278.394711096513});

    const double out = xsf::struve_l(v, x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-12;
    CAPTURE(v, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
