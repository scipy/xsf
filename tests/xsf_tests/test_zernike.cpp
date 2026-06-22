#include "../testing_utils.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <tuple>

#include <xsf/zernike.h>

namespace {

void require_close(double value, double expected, double tol = 1e-14) {
    CAPTURE(value, expected, tol);
    REQUIRE(std::abs(value - expected) <= tol);
}

double factorial(std::ptrdiff_t n) {
    double value = 1.0;

    for (std::ptrdiff_t i = 2; i <= n; ++i) {
        value *= static_cast<double>(i);
    }

    return value;
}

double zernike_radial_factorial_reference(std::ptrdiff_t n, std::ptrdiff_t m, double rho) {
    m = m < 0 ? -m : m;
    double value = 0.0;
    double sign = 1.0;

    for (std::ptrdiff_t k = 0; k <= (n - m) / 2; ++k) {
        const double coeff =
            sign * factorial(n - k) /
            (factorial(k) * factorial((n + m) / 2 - k) * factorial((n - m) / 2 - k));

        value += coeff * std::pow(rho, static_cast<double>(n - 2 * k));
        sign = -sign;
    }

    return value;
}

} // namespace

TEST_CASE("zernike radial low-order exact polynomials", "[zernike][xsf_tests]") {
    const double rho = 0.25;
    const double rho2 = rho * rho;
    const double rho4 = rho2 * rho2;

    require_close(xsf::eval_zernike_radial(0, 0, rho), 1.0);
    require_close(xsf::eval_zernike_radial(1, 1, rho), rho);
    require_close(xsf::eval_zernike_radial(2, 0, rho), 2.0 * rho2 - 1.0);
    require_close(xsf::eval_zernike_radial(2, 2, rho), rho2);
    require_close(xsf::eval_zernike_radial(4, 0, rho), 6.0 * rho4 - 6.0 * rho2 + 1.0);
    require_close(xsf::eval_zernike_radial(4, 2, rho), 4.0 * rho4 - 3.0 * rho2);
    require_close(xsf::eval_zernike_radial(4, 4, rho), rho4);
}

TEST_CASE("zernike radial identities", "[zernike][xsf_tests]") {
    for (std::ptrdiff_t n = 0; n <= 20; ++n) {
        require_close(xsf::eval_zernike_radial(n, n, 0.3), std::pow(0.3, static_cast<double>(n)));
        require_close(xsf::eval_zernike_radial(n, n, 1.0), 1.0);

        for (std::ptrdiff_t m = 0; m <= n; ++m) {
            if (((n - m) & 1) != 0) {
                continue;
            }

            require_close(xsf::eval_zernike_radial(n, m, 1.0), 1.0, 2e-13);
            require_close(xsf::eval_zernike_radial(n, -m, 0.6), xsf::eval_zernike_radial(n, m, 0.6));
        }
    }
}

TEST_CASE("zernike radial agrees with factorial definition", "[zernike][xsf_tests]") {
    using test_case = std::tuple<std::ptrdiff_t, std::ptrdiff_t, double>;
    auto [n, m, rho] = GENERATE(
        test_case{6, 0, 0.2}, test_case{6, 2, 0.4}, test_case{8, 4, 0.7}, test_case{10, 0, 0.3},
        test_case{12, 6, 0.9}, test_case{14, 2, 0.5}, test_case{18, 8, 0.75}
    );

    const double output = xsf::eval_zernike_radial(n, m, rho);
    const double expected = zernike_radial_factorial_reference(n, m, rho);
    const double error = xsf::extended_relative_error(output, expected);

    CAPTURE(n, m, rho, output, expected, error);
    REQUIRE(error <= 5e-13);
}

TEST_CASE("zernike radial barmak and kintner methods agree", "[zernike][xsf_tests]") {
    using test_case = std::tuple<std::ptrdiff_t, std::ptrdiff_t, double>;
    auto [n, m, rho] = GENERATE(
        test_case{0, 0, 0.0}, test_case{2, 0, 0.25}, test_case{5, 1, 0.5}, test_case{20, 0, 0.7},
        test_case{32, 12, 0.8}, test_case{50, 10, 0.7}, test_case{80, 30, 0.6}, test_case{100, 50, 0.9}
    );

    const double barmak = xsf::eval_zernike_radial_barmak(n, m, rho);
    const double kintner = xsf::eval_zernike_radial_kintner(n, m, rho);
    const double error = xsf::extended_relative_error(barmak, kintner);

    CAPTURE(n, m, rho, barmak, kintner, error);
    REQUIRE(error <= 5e-12);
}

TEST_CASE("zernike radial domain handling", "[zernike][xsf_tests]") {
    REQUIRE(std::isnan(xsf::eval_zernike_radial(-1, 0, 0.5)));
    REQUIRE(std::isnan(xsf::eval_zernike_radial(2, 3, 0.5)));
    REQUIRE(std::isnan(xsf::eval_zernike_radial(2, 1, 0.5)));
    REQUIRE(std::isnan(xsf::eval_zernike_radial(2, 0, -0.1)));
    REQUIRE(std::isnan(xsf::eval_zernike_radial(2, 0, 1.1)));
    REQUIRE(std::isnan(xsf::eval_zernike_radial(0, 0, std::numeric_limits<double>::quiet_NaN())));
}
