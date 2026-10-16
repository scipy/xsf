/* Regression tests for xsf::betaln at large arguments (scipy/scipy#26423).
 *
 * Reference values are lgamma(a) + lgamma(b) - lgamma(a + b) evaluated with
 * mpmath at a precision where a + b is exact.
 */

#include "../testing_utils.h"

#include <xsf/beta.h>

TEST_CASE("betaln large arguments with ratio above the asymptotic threshold", "[betaln][xsf_tests]") {
    // Before the fix lbeta_asymp overflowed b*(1-b), b*b*(1-b)*(1-b) or a*a*a
    // and returned NaN (inf/inf) for every case in the first group.
    constexpr double rtol = 1e-13;
    struct Case {
        double a, b, expected;
    };
    const std::vector<Case> cases = {
        {1e100, 1e107, -1.711809570095832e+101},
        {1e150, 1e157, -1.711809570095832e+151},
        {1e200, 1e230, -7.007755278982136e+201},
        {1e7, 1e100, -2151404143.624572},
        // ratio below the threshold: the non-asymptotic path, unchanged by the fix
        {1e10, 1e20, -240258509310.49857},
        {1e20, 1e30, -2.4025850929990457e+21},
        {2.0, 3.0, -2.4849066497880004},
    };
    for (const auto &c : cases) {
        const auto output = xsf::betaln(c.a, c.b);
        const auto rel_error = xsf::extended_relative_error(output, c.expected);
        CAPTURE(c.a, c.b, output, c.expected, rtol, rel_error);
        REQUIRE(std::isfinite(output));
        REQUIRE(rel_error <= rtol);
    }
}
