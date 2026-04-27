#include "../testing_utils.h"
#include <xsf/stats.h>

TEST_CASE("bivariate normal CDF test", "[bivariate_normal_cdf][xsf_tests]") {

    SECTION("bivariate normal CDF infinite inputs") {
        // test cases with infinite inputs
        using test_case = std::tuple<double, double, double, double, double>;
        auto [dh, dk, r, expected, rtol] = GENERATE(
            test_case{std::numeric_limits<double>::infinity(), 0.0, 0.5, 0.0, 1e-13},
            test_case{
                -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), 0.5, 1.0, 1e-13
            },
            test_case{-std::numeric_limits<double>::infinity(), 1.0, 0.5, xsf::ndtr(-1.0), 1e-13},
            test_case{1.0, -std::numeric_limits<double>::infinity(), 0.5, xsf::ndtr(-1.0), 1e-13}
        );
        const double output = xsf::bivariate_normal_cdf(dh, dk, r);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(dh, dk, r, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("bivariate normal CDF analytical dh=0, dk=0") {
        // bivariate_normal_cdf(0, 0, r) = 0.25 + asin(r)/(2*pi)
        double dh = 0.0;
        double dk = 0.0;
        double r = GENERATE(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0);
        double expected = 0.25 + std::asin(r) / (2 * M_PI);
        const double output = xsf::bivariate_normal_cdf(dh, dk, r);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        double rtol = 1e-13;
        CAPTURE(dh, dk, r, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("bivariate normal CDF r=0 independence") {
        // bivariate_normal_cdf(h, k, 0) = ndtr(-h) * ndtr(-dk)
        double dh = GENERATE(-1.0, 0.0, 1.0, 2.0);
        double dk = GENERATE(-1.0, 0.0, 1.0, 2.0);
        double r = 0.0;
        const double output = xsf::bivariate_normal_cdf(dh, dk, r);
        const double expected = xsf::ndtr(-dh) * xsf::ndtr(-dk);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        double rtol = 1e-13;
        CAPTURE(dh, dk, r, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("bivariate normal CDF symmetry h,k") {
        // bivariate_normal_cdf(h, k, r) == bivariate_normal_cdf(k, h, r)
        REQUIRE(xsf::bivariate_normal_cdf(1.0, 2.0, 0.5) == xsf::bivariate_normal_cdf(2.0, 1.0, 0.5));
    }

    SECTION("bivariate normal CDF complementarity") {
        // bivariate_normal_cdf(h, k, r) + bivariate_normal_cdf(h, -k, -r) == ndtr(-h)
        double dh = GENERATE(-1.0, 0.0, 1.0, 2.0);
        double dk = GENERATE(-1.0, 0.0, 1.0, 2.0);
        double r = GENERATE(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0);
        double output = xsf::bivariate_normal_cdf(dh, dk, r) + xsf::bivariate_normal_cdf(dh, -dk, -r);
        double expected = xsf::ndtr(-dh);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        double rtol = 1e-13;
        CAPTURE(dh, dk, r, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
