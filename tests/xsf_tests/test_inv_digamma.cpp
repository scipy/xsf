#include "../testing_utils.h"
#include <xsf/cephes/psi.h>
#include <xsf/inv_digamma.h>

#include <cmath>
#include <limits>

TEST_CASE("inv_digamma round-trip", "[inv_digamma][xsf_tests]") {
    const double rtol = 100 * std::numeric_limits<double>::epsilon();

    const std::vector<double> xs = {1e-10, 1e-5, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0};
    for (double x : xs) {
        const double y = xsf::cephes::psi(x);
        const auto result = xsf::inv_digamma(y);
        const auto relative_error = xsf::extended_relative_error(result, x);
        CAPTURE(x, y, result, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}

TEST_CASE("inv_digamma special values", "[inv_digamma][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    REQUIRE(std::isnan(xsf::inv_digamma(nan)));
    REQUIRE(xsf::inv_digamma(inf) == inf);
    REQUIRE(xsf::inv_digamma(-inf) == 0.0);
}
