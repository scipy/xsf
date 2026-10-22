#include "../testing_utils.h"
#include <xsf/beta.h>

TEST_CASE("betaln large arguments precision and stability", "[betaln][xsf_tests]") {
    constexpr double rtol = 1e-13;
    struct Case {
        double a, b, expected;
    };
    const std::vector<Case> cases = {
        // Cases that formerly returned NaN due to intermediate overflow
        {1e100, 1e107, -1.711809570095832e+101},
        {1e150, 1e157, -1.711809570095832e+151},
        {1e200, 1e230, -7.007755278982136e+201},
        {1e7, 1e100, -2151404143.624572},
        // Cases that formerly suffered from catastrophic cancellation (~1e-9 error)
        {2.2350235285263674e+123, 3.029877315654731e+117, -4.396727951061843e+118},
        {1.724473080349204e+106, 7.006411699535731e+100, -9.39812631468069e+101},
        {1e10, 1e20, -240258509310.49857},
        {1e20, 1e30, -2.4025850929990457e+21},
        // Standard range sanity check
        {2.0, 3.0, -2.4849066497880004},
        {1e5, 1e5, -1.386339270613481e+05},
    };

    for (const auto &c : cases) {
        const auto output = xsf::betaln(c.a, c.b);
        const auto rel_error = xsf::extended_relative_error(output, c.expected);
        CAPTURE(c.a, c.b, output, c.expected, rtol, rel_error);
        REQUIRE(std::isfinite(output));
        REQUIRE(rel_error <= rtol);
    }
}
