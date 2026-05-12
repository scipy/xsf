#include "../testing_utils.h"
#include <xsf/stats.h>
#include <cmath>
#include <limits>

TEST_CASE("gdtria basic", "[gdtria][xsf_tests]") {
    // From scipy docs: gdtr(1.2, 3.4, 5.6) = 0.9437808744201949
    // Verify inverse: gdtria(p, 3.4, 5.6) = 1.2
    //Run for both float and double precision
    const double p = 0.9437808744201949;
    const double a = 1.2;
    const double b = 3.4;
    const double x = 5.6;

    {
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto result = xsf::gdtria(p, b, x);
        const auto relative_error = xsf::extended_relative_error(result, a);
        CAPTURE(p, b, x, result, a, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }

    {
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto result = xsf::gdtria(static_cast<float>(p), static_cast<float>(b), static_cast<float>(x));
        const auto relative_error = static_cast<float>(xsf::extended_relative_error(result, static_cast<float>(a)));
        CAPTURE(p, b, x, result, a, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}

TEST_CASE("gdtria parameter boundaries", "[gdtria][xsf_tests]") {
    REQUIRE(std::isnan(xsf::gdtria(0.5, 2.0, 0.0)));
    REQUIRE(std::isnan(xsf::gdtria(0.5, -1.0, 5.0)));
    REQUIRE(std::isnan(xsf::gdtria(0.5, 2.0, -1.0)));
}

TEST_CASE("gdtrix basic", "[gdtrix][xsf_tests]") {
    // From scipy docs: gdtr(1.2, 3.4, 5.6) = 0.9437808744201949
    // Verify inverse: gdtrix(1.2, 3.4, p) = 5.6
    //Run for both float and double precision
    const double a = 1.2;
    const double b = 3.4;
    const double p = 0.9437808744201949;
    const double x = 5.6;

    {
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto result = xsf::gdtrix(a, b, p);
        const auto relative_error = xsf::extended_relative_error(result, x);
        CAPTURE(a, b, p, result, x, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }

    {
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto result = xsf::gdtrix(static_cast<float>(a), static_cast<float>(b), static_cast<float>(p));
        const auto relative_error = static_cast<float>(xsf::extended_relative_error(result, static_cast<float>(x)));
        CAPTURE(a, b, p, result, x, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}

TEST_CASE("gdtrix probability boundaries", "[gdtrix][xsf_tests]") {
    REQUIRE(std::isnan(xsf::gdtrix(0.0, 0.0, 0.5)));
    REQUIRE(xsf::gdtrix(1.0, 1.0, 0.0) == 0.0);
    REQUIRE(std::isinf(xsf::gdtrix(1.0, 1.0, 1.0)));
    REQUIRE(std::isnan(xsf::gdtrix(1.0, 1.0, -0.1)));
    REQUIRE(std::isnan(xsf::gdtrix(1.0, 1.0, 1.1)));
}