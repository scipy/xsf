#include "../testing_utils.h"
#include <xsf/cephes/psi.h>
#include <xsf/digamma_inv.h>

#include <cmath>
#include <limits>

TEST_CASE("digamma_inv round-trip", "[digamma_inv][xsf_tests]") {
    const double rtol = 100 * std::numeric_limits<double>::epsilon();

    const std::vector<double> xs = {1e-10, 1e-5, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0};
    for (double x : xs) {
        const double y = xsf::cephes::psi(x);
        const auto result = xsf::digamma_inv(y);
        const auto relative_error = xsf::extended_relative_error(result, x);
        CAPTURE(x, y, result, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}

TEST_CASE("digamma_inv special values", "[digamma_inv][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double log_max_double = std::log(std::numeric_limits<double>::max());

    REQUIRE(std::isnan(xsf::digamma_inv(nan)));
    REQUIRE(xsf::digamma_inv(inf) == inf);
    REQUIRE(xsf::digamma_inv(-inf) == 0.0);
    REQUIRE(xsf::digamma_inv(std::nextafter(log_max_double, inf)) == inf);
}
