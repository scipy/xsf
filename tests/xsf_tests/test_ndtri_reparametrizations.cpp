#include "../testing_utils.h"
#include <xsf/ndtri_reparametrizations.h>

#include <cmath>
#include <limits>

TEST_CASE("nrdtrimn basic", "[nrdtrimn][xsf_tests]") {
    // Test from the scipy docs example: mean=3, std=2, x=6 => p=ndtr((6-3)/2)
    // Run for both double and float precision.
    const double mean = 3.0;
    const double std = 2.0;
    const double x = 6.0;
    const double p = 0.9331927987311419;

    {
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto result = xsf::nrdtrimn(p, std, x);
        const auto relative_error = xsf::extended_relative_error(result, mean);
        CAPTURE(p, std, x, result, mean, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }

    {
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto result = xsf::nrdtrimn(static_cast<float>(p), static_cast<float>(std), static_cast<float>(x));
        const auto relative_error = static_cast<float>(xsf::extended_relative_error(result, static_cast<float>(mean)));
        CAPTURE(p, std, x, result, mean, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}

TEST_CASE("nrdtrimn NaN inputs", "[nrdtrimn][xsf_tests]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    REQUIRE(std::isnan(xsf::nrdtrimn(nan, 1.0, 0.0)));  // NaN p
    REQUIRE(std::isnan(xsf::nrdtrimn(0.5, nan, 0.0)));  // NaN std
    REQUIRE(std::isnan(xsf::nrdtrimn(0.5, 1.0, nan)));  // NaN x
    REQUIRE(std::isnan(xsf::nrdtrimn(0.5, 0.0, 0.0)));  // std == 0
    REQUIRE(std::isnan(xsf::nrdtrimn(0.5, -1.0, 0.0))); // std < 0
    REQUIRE(std::isnan(xsf::nrdtrimn(0.0, 1.0, 0.0)));  // p == 0
    REQUIRE(std::isnan(xsf::nrdtrimn(1.0, 1.0, 0.0)));  // p == 1
}

TEST_CASE("nrdtrisd basic", "[nrdtrisd][xsf_tests]") {
    // Test from the scipy docs example: mean=3, std=2, x=6 => p=ndtr((6-3)/2)
    // Run for both double and float precision.
    const double mean = 3.0;
    const double std = 2.0;
    const double x = 6.0;
    const double p = 0.9331927987311419;

    {
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto result = xsf::nrdtrisd(mean, p, x);
        const auto relative_error = xsf::extended_relative_error(result, std);
        CAPTURE(mean, p, x, result, std, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }

    {
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto result = xsf::nrdtrisd(static_cast<float>(mean), static_cast<float>(p), static_cast<float>(x));
        const auto relative_error = static_cast<float>(xsf::extended_relative_error(result, static_cast<float>(std)));
        CAPTURE(mean, p, x, result, std, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}

TEST_CASE("nrdtrisd NaN inputs", "[nrdtrisd][xsf_tests]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    REQUIRE(std::isnan(xsf::nrdtrisd(nan, 0.5, 0.0)));  // NaN mean
    REQUIRE(std::isnan(xsf::nrdtrisd(0.0, nan, 0.0)));  // NaN p
    REQUIRE(std::isnan(xsf::nrdtrisd(0.0, 0.5, nan)));  // NaN x
    REQUIRE(std::isnan(xsf::nrdtrisd(0.0, 0.0, 1.0)));  // p == 0
    REQUIRE(std::isnan(xsf::nrdtrisd(0.0, 1.0, 1.0)));  // p == 1
}