#include "../testing_utils.h"

#include <xsf/wright.h>

TEST_CASE("wrightomega nan complex", "[wrightomega][xsf_tests]") {
    auto p = GENERATE(
        std::complex<double>{NAN, 0.0}, std::complex<double>{0.0, NAN}, std::complex<double>{NAN, NAN},
        std::complex<double>{NAN, 1.0}, std::complex<double>{1.0, NAN}
    );

    const auto res = xsf::wrightomega(p);
    CAPTURE(p, res);
    REQUIRE(std::isnan(res.real()));
    REQUIRE(std::isnan(res.imag()));
}

TEST_CASE("wrightomega nan real", "[wrightomega][xsf_tests]") {
    auto res = xsf::wrightomega(NAN);
    CAPTURE(res);
    REQUIRE(std::isnan(res));
}

TEST_CASE("wrightomega inf branch", "[wrightomega][xsf_tests]") {
    using test_case = std::tuple<std::complex<double>, std::complex<double>>;
    auto [p, expected] = GENERATE(
        test_case{std::complex<double>{-INFINITY, M_PI / 4}, std::complex<double>{0.0, 0.0}},
        test_case{std::complex<double>{-INFINITY, -M_PI / 4}, std::complex<double>{0.0, -0.0}},
        test_case{std::complex<double>{-INFINITY, 3 * M_PI / 4}, std::complex<double>{-0.0, 0.0}},
        test_case{std::complex<double>{-INFINITY, -3 * M_PI / 4}, std::complex<double>{-0.0, -0.0}}
    );

    const auto res = xsf::wrightomega(p);
    // Compare real and imaginary parts separately to preserve signed zero checks.
    CAPTURE(p, res, expected);
    REQUIRE(res.real() == expected.real());
    REQUIRE(res.imag() == expected.imag());
}

TEST_CASE("wrightomega inf", "[wrightomega][xsf_tests]") {
    auto p = GENERATE(
        std::complex<double>{INFINITY, 10.0}, std::complex<double>{-INFINITY, 10.0},
        std::complex<double>{10.0, INFINITY}, std::complex<double>{10.0, -INFINITY}
    );

    const auto res = xsf::wrightomega(p);
    CAPTURE(p, res);
    REQUIRE(res == p);
}

TEST_CASE("wrightomega singular", "[wrightomega][xsf_tests]") {
    auto p = GENERATE(std::complex<double>{-1.0, M_PI}, std::complex<double>{-1.0, -M_PI});

    const auto res = xsf::wrightomega(p);
    CAPTURE(p, res);
    REQUIRE(res.real() == -1.0);
    REQUIRE(res.imag() == 0.0);
    REQUIRE(!std::signbit(res.imag()));
}

TEST_CASE("wrightomega real infinities", "[wrightomega][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(test_case{-INFINITY, 0.0}, test_case{INFINITY, INFINITY});

    const auto output = xsf::wrightomega(x);
    CAPTURE(x, output, desired);
    REQUIRE(output == desired);
}

TEST_CASE("wrightomega real series crossover", "[wrightomega][xsf_tests]") {
    constexpr double desired_error = 2 * std::numeric_limits<double>::epsilon();
    constexpr double crossover = 1e20;

    using test_case = std::tuple<double, double>;
    auto [x, expected] = GENERATE(
        test_case{std::nextafter(crossover, -INFINITY), 99999999999999983569.948},
        test_case{std::nextafter(crossover, INFINITY), 100000000000000016337.948}
    );

    const double output = xsf::wrightomega(x);
    const double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(x, output, expected, desired_error, rel_error);
    REQUIRE(rel_error <= desired_error);
}

TEST_CASE("wrightomega exp approximation crossover", "[wrightomega][xsf_tests]") {
    constexpr double desired_error = 2 * std::numeric_limits<double>::epsilon();
    constexpr double crossover = -50.0;

    using test_case = std::tuple<double, double>;
    auto [x, expected] = GENERATE(
        test_case{std::nextafter(crossover, INFINITY), 1.9287498479639314876e-22},
        test_case{std::nextafter(crossover, -INFINITY), 1.9287498479639040784e-22}
    );

    const double output = xsf::wrightomega(x);
    const double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(x, output, expected, desired_error, rel_error);
    REQUIRE(rel_error <= desired_error);
}

TEST_CASE("wrightomega real versus complex", "[wrightomega][xsf_tests]") {
    constexpr double rtol = 1e-14;
    const auto xs = linspace(-500.0, 500.0, 1001);

    for (double x : xs) {
        const double output = xsf::wrightomega(x);
        const double expected = std::real(xsf::wrightomega(std::complex<double>{x, 0.0}));
        const double rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
