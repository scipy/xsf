#include "../testing_utils.h"
#include "xsf/config.h"
#include <xsf/iv_ratio.h>

TEST_CASE("iv_ratioinv round-trip single-double", "[iv_ratioinv][xsf_tests]") {
    const std::vector<double> xs = logspace<double>(-10, 7, 100);
    const std::vector<double> vs = logspace<double>(0, 6, 100);
    for (double x : xs) {
        for (double v : vs) {
            const double y_forward = xsf::iv_ratio(v, x);
            const double y = (y_forward < 0.5) ? y_forward : (1.0 - xsf::iv_ratio_c(v, x));
            if ((y == 0.) || (y == 1.)) {
                continue;
            }
            const auto result = xsf::iv_ratioinv(v, y);
            const auto relative_error = xsf::extended_relative_error(result, x);
            double rtol = 1e-9;
            CAPTURE(v, x, y, result, rtol, relative_error);
            REQUIRE(relative_error <= rtol);

            const float v_f = static_cast<float>(v);
            const float x_f = static_cast<float>(x);
            const float y_forward_f = xsf::iv_ratio(v_f, x_f);
            const float y_f = (y_forward_f < 0.5f) ? y_forward_f : (1.0f - xsf::iv_ratio_c(v_f, x_f));
            if ((y_f < 0.0001f) || (y_f > 0.9999f)) {
                // Skip the tails in single precision because precision is not sufficient
                continue;
            }
            const auto result_f = xsf::iv_ratioinv(v_f, y_f);
            const auto relative_error_f = xsf::extended_relative_error(result_f, x_f);
            float rtol_f = 5e-4f;
            CAPTURE(v, x, y_f, result_f, rtol_f, relative_error_f);
            REQUIRE(relative_error_f <= rtol_f);
        }
    }
}

TEST_CASE("iv_ratioinv arbitrary precision", "[iv_ratioinv][xsf_tests]") {
    // Test with arbitrary precision values for v, x, and iv_ratio(v,x)
    // Reference values for iv_ratio were computed with the Python library mpmath.
    // from mpmath import mp
    // mp.dps = 1000
    // def iv_ratio(v, x):
    //     v = mp.mpf(v)
    //     x = mp.mpf(x)
    //     result = mp.besseli(v, x) / mp.besseli(v - 1, x)
    //     return float(result)

    using test_case = std::tuple<double, double, double>;
    auto [v, x, iv_ratio_ref] = GENERATE(
        test_case{10, 100000, 0.9999050040375403}, test_case{1, 0.1, 0.04993760398793892},
        test_case{0.5, 0.001, 0.0009999996666668}, test_case{1000, 500, 0.23607911616885813},
        test_case{100, 1e-15, 5e-18}, test_case{1e5, 5, 2.4999999984375157e-05},
        test_case{1e5, 1e5, 0.41421399130458997}, test_case{1e5, 1e10, 0.999990000099999},
        test_case{0.5, 18.714973875118524, std::nextafter(1.0, 0.0)},
        test_case{0.5, std::numeric_limits<double>::min(), std::numeric_limits<double>::min()}
    );
    double x_result = xsf::iv_ratioinv(v, iv_ratio_ref);
    const auto rel_error = xsf::extended_relative_error(x_result, x);

    CAPTURE(v, x, iv_ratio_ref, x_result, rel_error);
    REQUIRE(rel_error <= 1e-10);
}

TEST_CASE("iv_ratioinv nan propagation", "[iv_ratioinv][xsf_tests]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(std::isnan(xsf::iv_ratioinv(nan, 0.5)));
    REQUIRE(std::isnan(xsf::iv_ratioinv(1.0, nan)));
}

TEST_CASE("iv_ratioinv domain error", "[iv_ratioinv][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    REQUIRE(std::isnan(xsf::iv_ratioinv(0.4, 0.5)));  // v < 0.5
    REQUIRE(std::isnan(xsf::iv_ratioinv(1.0, -1.0))); // y < 0
    REQUIRE(std::isnan(xsf::iv_ratioinv(1.0, 2.0)));  // y > 1
    REQUIRE(std::isnan(xsf::iv_ratioinv(inf, 0.5)));  // v = inf
}

TEST_CASE("iv_ratioinv domain boundary", "[iv_ratioinv][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    REQUIRE(xsf::iv_ratioinv(1.0, 0.0) == 0.0); // y = 0
    REQUIRE(xsf::iv_ratioinv(1.0, 1.0) == inf); // y = 1
}

TEST_CASE("iv_ratioinv v = 0.5 roundtrip", "[iv_ratioinv][xsf_tests]") {
    const std::vector<double> xs = logspace<double>(-50, 2, 500);
    double threshold = 1.0 - 1e-12; // threshold to avoid numerical issues near y=1
    for (double x : xs) {
        const double y_forward = xsf::iv_ratio(0.5, x);
        const double y = (y_forward < 0.5) ? y_forward : (1.0 - xsf::iv_ratio_c(0.5, x));
        if ((y == 0.) || (y > threshold)) {
            continue;
        }
        const double x_result = xsf::iv_ratioinv(0.5, y);
        const auto relative_error = xsf::extended_relative_error(x_result, x);
        double rtol = 1e-7;
        CAPTURE(x, y, x_result, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}
