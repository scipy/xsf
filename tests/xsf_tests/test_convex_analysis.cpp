#include "../testing_utils.h"

#include <xsf/convex_analysis.h>

namespace {

constexpr double rtol = 1e-13;
constexpr double abs_tol = 1e-13;

} // namespace

TEST_CASE("entr", "[entr][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4563
    const std::vector<double> values = {0, 0.5, 1.0, std::numeric_limits<double>::infinity()};
    const std::vector<double> signs = {-1, 1};

    for (double sign : signs) {
        for (double value : values) {
            double x = sign * value;
            double desired;
            if (x < 0) {
                desired = -std::numeric_limits<double>::infinity();
            } else if (x == 0) {
                desired = 0;
            } else {
                desired = -x * std::log(x);
            }

            auto out = xsf::entr(x);
            auto abs_error = xsf::extended_absolute_error(out, desired);
            auto rel_error = xsf::extended_relative_error(out, desired);
            CAPTURE(x, out, desired, abs_error, rel_error, abs_tol, rtol);
            REQUIRE((abs_error <= abs_tol || rel_error <= rtol));
        }
    }
}

TEST_CASE("kl_div", "[kl_div][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4579
    const std::vector<double> values = {0, 0.5, 1.0};
    const std::vector<double> signs = {-1, 1};

    for (double sgna : signs) {
        for (double va : values) {
            for (double sgnb : signs) {
                for (double vb : values) {
                    double x = sgna * va;
                    double y = sgnb * vb;
                    double desired;
                    if (x < 0 || y < 0 || (y == 0 && x != 0)) {
                        desired = std::numeric_limits<double>::infinity();
                    } else if ((std::isinf(x) && x > 0) || (std::isinf(y) && y > 0)) {
                        desired = std::numeric_limits<double>::infinity();
                    } else if (x == 0) {
                        desired = y;
                    } else {
                        desired = x * std::log(x / y) - x + y;
                    }
                    auto out = xsf::kl_div(x, y);
                    auto abs_error = xsf::extended_absolute_error(out, desired);
                    auto rel_error = xsf::extended_relative_error(out, desired);
                    CAPTURE(x, y, out, desired, abs_error, rel_error, abs_tol, rtol);
                    REQUIRE((abs_error <= abs_tol || rel_error <= rtol));
                }
            }
        }
    }
}

TEST_CASE("rel_entr", "[rel_entr][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4601
    const std::vector<double> values = {0, 0.5, 1.0};
    const std::vector<double> signs = {-1, 1};

    for (double sgna : signs) {
        for (double va : values) {
            for (double sgnb : signs) {
                for (double vb : values) {
                    double x = sgna * va;
                    double y = sgnb * vb;
                    double desired;
                    if (x > 0 && y > 0) {
                        desired = x * std::log(x / y);
                    } else if (x == 0 && y >= 0) {
                        desired = 0;
                    } else {
                        desired = std::numeric_limits<double>::infinity();
                    }
                    auto out = xsf::rel_entr(x, y);
                    auto abs_error = xsf::extended_absolute_error(out, desired);
                    auto rel_error = xsf::extended_relative_error(out, desired);
                    CAPTURE(x, y, out, desired, abs_error, rel_error, abs_tol, rtol);
                    REQUIRE((abs_error <= abs_tol || rel_error <= rtol));
                }
            }
        }
    }
}

TEST_CASE("rel_entr gh-20710 near zero", "[rel_entr][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4619
    const std::vector<std::tuple<double, double, double>> cases = {
        {0.9456657713430001, 0.9456657713430094, -9.325873406851269e-15},
        {0.48066098564791515, 0.48066098564794774, -3.258504577274724e-14},
        {0.786048657854401, 0.7860486578542367, 1.6431300764454033e-13},
    };

    for (auto [x, y, desired] : cases) {
        auto out = xsf::rel_entr(x, y);
        auto rel_error = xsf::extended_relative_error(out, desired);
        CAPTURE(x, y, out, desired, rel_error, rtol);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("rel_entr gh-20710 overflow", "[rel_entr][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4638
    const std::vector<std::tuple<double, double, double>> cases = {
        {4, 2.22e-308, 2839.139983229607},
        {1e-200, 1e+200, -9.210340371976183e-198},
        {2.22e-308, 1e15, -1.6493212008074475e-305},
    };

    for (auto [x, y, desired] : cases) {
        auto out = xsf::rel_entr(x, y);
        auto rel_error = xsf::extended_relative_error(out, desired);
        CAPTURE(x, y, out, desired, rel_error, rtol);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("huber", "[huber][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4660
    REQUIRE(xsf::huber(-1, 1.5) == std::numeric_limits<double>::infinity());
    REQUIRE(xsf::extended_absolute_error(xsf::huber(2, 1.5), 0.5 * 1.5 * 1.5) <= abs_tol);
    REQUIRE(xsf::extended_absolute_error(xsf::huber(2, 2.5), 2 * (2.5 - 0.5 * 2)) <= abs_tol);

    const std::vector<std::tuple<double, double>> cases = {
        {-1.25, 0.5}, {0, 1}, {0.25, -0.1}, {0.25, -0.25}, {0.25, -0.5}, {1, 0.5}, {1, 1}, {1, 2}, {3, -4}, {5, 4},
    };

    for (auto [delta, r] : cases) {
        double desired;
        if (delta < 0) {
            desired = std::numeric_limits<double>::infinity();
        } else if (std::abs(r) < delta) {
            desired = 0.5 * r * r;
        } else {
            desired = delta * (std::abs(r) - 0.5 * delta);
        }
        auto out = xsf::huber(delta, r);
        auto abs_error = xsf::extended_absolute_error(out, desired);
        auto rel_error = xsf::extended_relative_error(out, desired);
        CAPTURE(delta, r, out, desired, abs_error, rel_error, abs_tol, rtol);
        REQUIRE((abs_error <= abs_tol || rel_error <= rtol));
    }
}

TEST_CASE("pseudo_huber", "[pseudo_huber][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4678
    auto deltas = linspace(-5.0, 5.0, 11);
    auto residuals = linspace(-5.0, 5.0, 11);
    std::vector<std::tuple<double, double>> cases(deltas.size());
    for (std::size_t i = 0; i < deltas.size(); ++i) {
        cases[i] = {deltas[i], residuals[i]};
    }
    cases.push_back({0, 0.5});
    cases.push_back({0.5, 0});

    for (auto [delta, r] : cases) {
        double desired;
        if (delta < 0) {
            desired = std::numeric_limits<double>::infinity();
        } else if (!delta || !r) {
            desired = 0;
        } else {
            desired = delta * delta * (std::sqrt(1 + (r / delta) * (r / delta)) - 1);
        }
        auto out = xsf::pseudo_huber(delta, r);
        auto abs_error = xsf::extended_absolute_error(out, desired);
        auto rel_error = xsf::extended_relative_error(out, desired);
        CAPTURE(delta, r, out, desired, abs_error, rel_error, abs_tol, rtol);
        REQUIRE((abs_error <= abs_tol || rel_error <= rtol));
    }
}

TEST_CASE("pseudo_huber small r", "[pseudo_huber][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4692
    auto out = xsf::pseudo_huber(1.0, 1e-18);
    auto desired = 5.0000000000000005e-37;
    auto rel_error = xsf::extended_relative_error(out, desired);
    CAPTURE(out, desired, rel_error, rtol);
    REQUIRE(rel_error <= rtol);
}
