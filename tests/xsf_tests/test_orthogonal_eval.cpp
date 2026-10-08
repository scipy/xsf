#include "../../include/xsf/config.h"
#include "../testing_utils.h"

#include <xsf/cpu/orthogonal_eval.h>
#include <xsf/evalpoly.h>
#include <xsf/gamma.h>
#include <xsf/orthogonal_eval.h>

#include <array>
#include <vector>

namespace {

double binom_int(double a, int k) {
    double out = 1.0;
    for (int j = 0; j < k; ++j) {
        out *= (a - j) / (j + 1.0);
    }
    return out;
}

std::vector<double> laguerre_coefficients(int n, double alpha) {
    std::vector<double> out(n + 1, 0.0);
    double scale = 1.0;
    for (int k = 0; k <= n; ++k) {
        out[k] = scale * binom_int(n + alpha, n - k);
        scale /= -(k + 1.0);
    }
    return out;
}

std::vector<double> linear_power(double c0, double c1, int n) {
    std::vector<double> out(n + 1, 0.0);
    for (int j = 0; j <= n; ++j) {
        out[j] = binom_int(n, j) * std::pow(c0, n - j) * std::pow(c1, j);
    }
    return out;
}

std::vector<double> multiply(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> out(a.size() + b.size() - 1, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            out[i + j] += a[i] * b[j];
        }
    }
    return out;
}

template <typename T>
T polyval(const std::vector<double> &coeffs, T x) {
    if (coeffs.size() == 1) {
        return coeffs[0];
    }

    const std::vector<double> reversed(coeffs.rbegin(), coeffs.rend());
    return xsf::evalpoly(reversed.data(), reversed.size() - 1, x);
}

double sample(double a, double b, int i) {
    constexpr double step = 0.7548776662466927;
    const double u = std::fmod(step * (i + 1), 1.0);
    return a + (b - a) * u;
}

template <typename CoefficientsFunc, typename EvalFunc>
void check_poly(
    CoefficientsFunc coefficients_func, EvalFunc eval_func, const std::vector<std::pair<double, double>> &param_ranges,
    std::pair<double, double> x_range, double rtol, int nn = 10, int nparam = 10, int nx = 10
) {
    for (int n = 0; n < nn; ++n) {
        const int ncases = param_ranges.empty() ? 1 : nparam;
        for (int ip = 0; ip < ncases; ++ip) {
            std::vector<double> params;
            params.reserve(param_ranges.size());
            for (std::size_t k = 0; k < param_ranges.size(); ++k) {
                params.push_back(sample(param_ranges[k].first, param_ranges[k].second, 17 * n + nparam * k + ip));
            }
            const auto coeffs = coefficients_func(n, params);

            for (int ix = 0; ix < nx; ++ix) {
                double x = sample(x_range.first, x_range.second, 31 * n + 13 * ip + ix);
                if (ix == 0) {
                    x = x_range.first;
                } else if (ix == 1) {
                    x = x_range.second;
                }

                const auto out = eval_func(n, params, x);
                const double expected = polyval(coeffs, x);
                const double error = xsf::extended_absolute_error(static_cast<double>(out), expected);
                const double tol = 1e-12 + rtol * std::abs(expected);
                CAPTURE(n, params, x, out, expected, error, tol);
                REQUIRE(error <= tol);
            }
        }
    }
}

template <typename CoefficientsFunc, typename EvalFunc>
void check_complex_poly(
    CoefficientsFunc coefficients_func, EvalFunc eval_func, const std::vector<std::pair<double, double>> &param_ranges,
    const std::vector<std::complex<double>> &xs, double rtol, int nn = 10, int nparam = 10
) {
    for (int n = 0; n < nn; ++n) {
        const int ncases = param_ranges.empty() ? 1 : nparam;
        for (int ip = 0; ip < ncases; ++ip) {
            std::vector<double> params;
            params.reserve(param_ranges.size());
            for (std::size_t k = 0; k < param_ranges.size(); ++k) {
                params.push_back(sample(param_ranges[k].first, param_ranges[k].second, 23 * n + nparam * k + ip));
            }
            const auto coeffs = coefficients_func(n, params);

            for (const auto &x : xs) {
                const auto out = eval_func(n, params, x);
                const auto expected = polyval(coeffs, x);
                const double error = xsf::extended_absolute_error(static_cast<std::complex<double>>(out), expected);
                const double tol = 1e-12 + rtol * std::abs(expected);
                CAPTURE(n, params, x, out, expected, error, tol);
                REQUIRE(error <= tol);
            }
        }
    }
}

template <typename IntDegreeEvalFunc, typename DoubleDegreeEvalFunc>
void check_recurrence(
    IntDegreeEvalFunc int_degree_eval, DoubleDegreeEvalFunc double_degree_eval,
    const std::vector<std::pair<double, double>> &param_ranges, std::pair<double, double> x_range, double rtol = 1e-8,
    int nn = 10, int nparam = 10, int nx = 10
) {
    for (int n = 0; n < nn; ++n) {
        const int ncases = param_ranges.empty() ? 1 : nparam;
        for (int ip = 0; ip < ncases; ++ip) {
            std::vector<double> params;
            params.reserve(param_ranges.size());
            for (std::size_t k = 0; k < param_ranges.size(); ++k) {
                params.push_back(sample(param_ranges[k].first, param_ranges[k].second, 19 * n + nparam * k + ip));
            }

            for (int ix = 0; ix < nx; ++ix) {
                double x = sample(x_range.first, x_range.second, 37 * n + 11 * ip + ix);
                if (ix == 0) {
                    x = x_range.first;
                } else if (ix == 1) {
                    x = x_range.second;
                }

                const double out = int_degree_eval(n, params, x);
                const double expected = double_degree_eval(n, params, x);
                const double error = xsf::extended_absolute_error(out, expected);
                const double tol = 1e-12 + rtol * std::abs(expected);
                CAPTURE(n, params, x, out, expected, error, tol);
                REQUIRE(error <= tol);
            }
        }
    }
}

} // namespace

TEST_CASE("eval_jacobi matches constructed polynomials", "[eval_jacobi][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L77-L80
    check_poly(
        [](int n, const std::vector<double> &params) {
            const double alpha = params[0];
            const double beta = params[1];
            std::vector<double> out(n + 1.0, 0.0);
            const double scale = std::ldexp(1.0, -n);
            for (int m = 0; m <= n; ++m) {
                const double c = scale * binom_int(n + alpha, m) * binom_int(n + beta, n - m);
                const auto term = multiply(linear_power(-1.0, 1.0, n - m), linear_power(1.0, 1.0, m));
                for (int j = 0; j <= n; ++j) {
                    out[j] += c * term[j];
                }
            }
            return out;
        },
        [](double n, const std::vector<double> &params, double x) {
            return xsf::eval_jacobi(n, params[0], params[1], x);
        },
        {{-0.99, 10.0}, {-0.99, 10.0}}, {-1.0, 1.0}, 1e-5
    );
}

TEST_CASE("eval_jacobi supports complex inputs", "[eval_jacobi][xsf_tests]") {
    // for complex inputs
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L77-L80
    check_complex_poly(
        [](double n, const std::vector<double> &params) {
            const double alpha = params[0];
            const double beta = params[1];
            std::vector<double> out(n + 1.0, 0.0);
            const double scale = std::ldexp(1.0, -n);
            for (int m = 0; m <= n; ++m) {
                const double c = scale * binom_int(n + alpha, m) * binom_int(n + beta, n - m);
                const auto term = multiply(linear_power(-1.0, 1.0, n - m), linear_power(1.0, 1.0, m));
                for (int j = 0; j <= n; ++j) {
                    out[j] += c * term[j];
                }
            }
            return out;
        },
        [](double n, const std::vector<double> &params, std::complex<double> x) {
            return xsf::eval_jacobi(n, params[0], params[1], x);
        },
        {{-0.99, 10.0}, {-0.99, 10.0}}, {{-0.75, -0.25}, {-0.25, 0.5}, {0.0, -0.5}, {0.5, 0.25}, {0.75, -0.75}}, 1e-5
    );
}

TEST_CASE("eval_jacobi matches SciPy complex<double> reference values", "[eval_jacobi][xsf_tests]") {
    using test_case = std::tuple<double, double, double, std::complex<double>, std::complex<double>>;
    auto [n, alpha, beta, x, expected] = GENERATE(
        test_case{2.0, 0.25, 1.5, {0.2, -0.3}, {-0.9910156249999997, 0.035625000000000136}},
        test_case{4.0, 2.5, -0.4, {-0.35, 0.6}, {5.171472426694353, 0.8954783753906292}},
        test_case{6.0, 0.75, 3.25, {0.9, -0.2}, {-7.36364426660156, 1.1136538789062485}}
    );

    const auto out = xsf::eval_jacobi(n, alpha, beta, x);
    const double error = xsf::extended_absolute_error(out, expected);
    const double tol = 1e-12 + 1e-12 * std::abs(expected);
    CAPTURE(n, alpha, beta, x, out, expected, error, tol);
    REQUIRE(error <= tol);
}

TEST_CASE("eval_jacobi matches SciPy complex<float> reference values", "[eval_jacobi][xsf_tests]") {
    using test_case = std::tuple<float, float, float, std::complex<float>, std::complex<float>>;
    auto [n, alpha, beta, x, expected] = GENERATE(
        test_case{2.0f, 0.25f, 1.5f, {0.2f, -0.3f}, {-0.99101567f, 0.035624996f}},
        test_case{4.0f, 2.5f, -0.4f, {-0.35f, 0.6f}, {5.171473f, 0.8954783f}},
        test_case{6.0f, 0.75f, 3.25f, {0.9f, -0.2f}, {-7.363643f, 1.1136553f}}
    );

    const auto out = xsf::eval_jacobi(n, alpha, beta, x);
    const double error = xsf::extended_absolute_error(out, expected);
    const double tol = 1e-6 + 1e-5 * std::abs(expected);
    CAPTURE(n, alpha, beta, x, out, expected, error, tol);
    REQUIRE(error <= tol);
}

TEST_CASE("eval_sh_jacobi matches constructed polynomials", "[eval_sh_jacobi][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L82-L85
    check_poly(
        [](int n, const std::vector<double> &params) {
            const double p = params[0];
            const double q = params[1];
            const double alpha = p - q;
            const double beta = q - 1.0;
            const double scale = 1.0 / xsf::binom(2.0 * n + p - 1.0, n);
            std::vector<double> out(n + 1, 0.0);
            for (int m = 0; m <= n; ++m) {
                const double c = scale * binom_int(n + alpha, m) * binom_int(n + beta, n - m);
                const auto term = multiply(linear_power(-1.0, 1.0, n - m), linear_power(0.0, 1.0, m));
                for (int j = 0; j <= n; ++j) {
                    out[j] += c * term[j];
                }
            }
            return out;
        },
        [](double n, const std::vector<double> &params, double x) {
            return xsf::eval_sh_jacobi(n, params[0], params[1], x);
        },
        {{1.0, 10.0}, {0.0, 1.0}}, {0.0, 1.0}, 1e-5
    );
}

TEST_CASE("eval_sh_jacobi for complex inputs", "[eval_sh_jacobi][xsf_tests]") {
    // for complex inputs
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L82-L85
    check_complex_poly(
        [](double n, const std::vector<double> &params) {
            const double p = params[0];
            const double q = params[1];
            const double alpha = p - q;
            const double beta = q - 1.0;
            const double scale = 1.0 / xsf::binom(2.0 * n + p - 1.0, n);
            std::vector<double> out(n + 1, 0.0);
            for (int m = 0; m <= n; ++m) {
                const double c = scale * binom_int(n + alpha, m) * binom_int(n + beta, n - m);
                const auto term = multiply(linear_power(-1.0, 1.0, n - m), linear_power(0.0, 1.0, m));
                for (int j = 0; j <= n; ++j) {
                    out[j] += c * term[j];
                }
            }
            return out;
        },
        [](double n, const std::vector<double> &params, std::complex<double> x) {
            return xsf::eval_sh_jacobi(n, params[0], params[1], x);
        },
        {{1.0, 10.0}, {0.0, 1.0}}, {{0.1, 0.2}, {0.25, -0.3}, {0.5, 0.4}, {0.75, -0.2}, {0.9, 0.1}}, 1e-5
    );
}

TEST_CASE("eval_sh_jacobi matches SciPy complex<double> reference values", "[eval_sh_jacobi][xsf_tests]") {
    using test_case = std::tuple<double, double, double, std::complex<double>, std::complex<double>>;
    auto [n, p, q, x, expected] = GENERATE(
        test_case{2.0, 1.25, 0.5, {0.2, 0.1}, {-0.05687782805429855, -0.030588235294117673}},
        test_case{4.0, 3.5, 0.25, {0.8, -0.15}, {-0.005872881652661048, -0.043102675807164974}},
        test_case{5.0, 2.0, 0.75, {-0.1, 0.4}, {0.09778388908617434, -0.07075225852272753}}
    );

    const auto out = xsf::eval_sh_jacobi(n, p, q, x);
    const double error = xsf::extended_absolute_error(out, expected);
    const double tol = 1e-12 + 1e-12 * std::abs(expected);
    CAPTURE(n, p, q, x, out, expected, error, tol);
    REQUIRE(error <= tol);
}

TEST_CASE("eval_sh_jacobi matches SciPy complex<float> reference values", "[eval_sh_jacobi][xsf_tests]") {
    using test_case = std::tuple<float, float, float, std::complex<float>, std::complex<float>>;
    auto [n, p, q, x, expected] = GENERATE(
        test_case{3.0f, 1.5f, 0.125f, {0.15f, 0.35f}, {0.074561186f, -0.05201661f}},
        test_case{5.0f, 4.25f, 0.625f, {0.6f, -0.25f}, {0.0034959577f, 0.009585135f}},
        test_case{6.0f, 2.75f, 0.9f, {-0.2f, 0.3f}, {-0.0819409f, 0.004271384f}}
    );

    const auto out = xsf::eval_sh_jacobi(n, p, q, x);
    const double error = xsf::extended_absolute_error(out, expected);
    const double tol = 1e-6 + 1e-5 * std::abs(expected);
    CAPTURE(n, p, q, x, out, expected, error, tol);
    REQUIRE(error <= tol);
}

TEST_CASE("eval_jacobi recurrence overload", "[eval_jacobi][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L185-L188
    check_recurrence(
        [](int n, const std::vector<double> &params, double x) { return xsf::eval_jacobi(n, params[0], params[1], x); },
        [](double n, const std::vector<double> &params, double x) {
            return xsf::eval_jacobi(n, params[0], params[1], x);
        },
        {{-0.99, 10.0}, {-0.99, 10.0}}, {-1.0, 1.0}
    );
}

TEST_CASE("eval_sh_jacobi recurrence overload", "[eval_sh_jacobi][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L190-L192
    check_recurrence(
        [](int n, const std::vector<double> &params, double x) {
            return xsf::eval_sh_jacobi(n, params[0], params[1], x);
        },
        [](double n, const std::vector<double> &params, double x) {
            return xsf::eval_sh_jacobi(n, params[0], params[1], x);
        },
        {{1.0, 10.0}, {0.0, 1.0}}, {0.0, 1.0}
    );
}

TEST_CASE("eval_jacobi alpha=-1 beta=1", "[eval_jacobi][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L293-L329
    using test_case = std::tuple<int, std::array<double, 11>>;
    // gh-7001 - expected values were computed with mathematica.
    auto [n, expected] = GENERATE(
        test_case{0, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}},
        test_case{1, {-2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0}},
        test_case{2, {3.0, 2.16, 1.44, 0.84, 0.36, 0.0, -0.24, -0.36, -0.36, -0.24, 0.0}},
        test_case{3, {-4.0, -1.98, -0.64, 0.14, 0.48, 0.5, 0.32, 0.06, -0.16, -0.22, 0.0}},
        test_case{4, {5.0, 1.332, -0.288, -0.658, -0.408, 0.0, 0.272, 0.282, 0.072, -0.148, 0.0}},
        test_case{5, {-6.0, -0.43308, 0.79104, 0.36876, -0.21312, -0.375, -0.14208, 0.15804, 0.19776, -0.04812, 0.0}}
    );

    for (std::size_t j = 0; j < expected.size(); ++j) {
        const double x = -1.0 + 0.2 * static_cast<double>(j);
        auto out = xsf::eval_jacobi(n, -1.0, 1.0, x);
        auto error = xsf::extended_absolute_error(out, expected[j]);
        auto tol = 1e-14 + 1e-10 * std::abs(expected[j]);
        CAPTURE(n, x, out, expected[j], error, tol);
        REQUIRE(error <= tol);

        out = xsf::eval_jacobi(static_cast<double>(n), -1.0, 1.0, x);
        error = xsf::extended_absolute_error(out, expected[j]);
        CAPTURE(n, x, out, expected[j], error, tol);
        REQUIRE(error <= tol);
    }
}

TEST_CASE("eval_jacobi alpha=-1 beta=-1", "[eval_jacobi][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L332-L383
    using test_case = std::tuple<int, std::array<double, 11>>;
    // gh-7001 - expected values were computed with mathematica.
    auto [n, expected] = GENERATE(
        test_case{0, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}},
        test_case{1, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        test_case{2, {0.0, -0.09, -0.16, -0.21, -0.24, -0.25, -0.24, -0.21, -0.16, -0.09, 0.0}},
        test_case{3, {0.0, 0.144, 0.192, 0.168, 0.096, 0.0, -0.096, -0.168, -0.192, -0.144, 0.0}},
        test_case{4, {0.0, -0.1485, -0.096, 0.0315, 0.144, 0.1875, 0.144, 0.0315, -0.096, -0.1485, 0.0}},
        test_case{5, {0.0, 0.10656, -0.04608, -0.15792, -0.13056, 0.0, 0.13056, 0.15792, 0.04608, -0.10656, 0.0}}
    );

    for (std::size_t j = 0; j < expected.size(); ++j) {
        const double x = -1.0 + 0.2 * static_cast<double>(j);
        auto out = xsf::eval_jacobi(n, -1.0, -1.0, x);
        auto error = xsf::extended_absolute_error(out, expected[j]);
        auto tol = 1e-14 + 1e-10 * std::abs(expected[j]);
        CAPTURE(n, x, out, expected[j], error, tol);
        REQUIRE(error <= tol);

        out = xsf::eval_jacobi(static_cast<double>(n), -1.0, -1.0, x);
        error = xsf::extended_absolute_error(out, expected[j]);
        CAPTURE(n, x, out, expected[j], error, tol);
        REQUIRE(error <= tol);
    }
}

TEST_CASE("eval_hermitenorm matches constructed polynomials", "[eval_hermitenorm][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_orthogonal_eval.py#L138-L140
    check_poly(
        [](int n, const std::vector<double> &) {
            std::vector<double> out(n + 1, 0.0);
            for (int m = 0; m <= n / 2; ++m) {
                const int power = n - 2 * m;
                out[power] = std::pow(-1.0, m) * xsf::gamma(n + 1.0) /
                             (std::pow(2.0, m) * xsf::gamma(m + 1.0) * xsf::gamma(power + 1.0));
            }
            return out;
        },
        [](int n, const std::vector<double> &, double x) { return xsf::eval_hermitenorm(n, x); }, {}, {-100.0, 100.0},
        1e-12
    );
}

TEST_CASE("eval_hermite matches constructed polynomials", "[eval_hermite][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_orthogonal_eval.py#L134-L136
    check_poly(
        [](int n, const std::vector<double> &) {
            std::vector<double> out(n + 1, 0.0);
            for (int m = 0; m <= n / 2; ++m) {
                const int power = n - 2 * m;
                out[power] = std::pow(-1.0, m) * xsf::gamma(n + 1.0) * std::pow(2.0, power) /
                             (xsf::gamma(m + 1.0) * xsf::gamma(power + 1.0));
            }
            return out;
        },
        [](int n, const std::vector<double> &, double x) { return xsf::eval_hermite(n, x); }, {}, {-100.0, 100.0}, 1e-12
    );
}

TEST_CASE("Hermite evaluators handle domain and NaN inputs", "[eval_hermite][eval_hermitenorm][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_orthogonal_eval.py#L244-L255
    REQUIRE(std::isnan(xsf::eval_hermite(-1, 1.0)));
    REQUIRE(std::isnan(xsf::eval_hermitenorm(-1, 1.0)));

    for (int n = 0; n <= 2; ++n) {
        for (double x : {0.0, 1.0, std::numeric_limits<double>::quiet_NaN()}) {
            CAPTURE(n, x);
            REQUIRE(std::isnan(xsf::eval_hermite(n, x)) == std::isnan(x));
            REQUIRE(std::isnan(xsf::eval_hermitenorm(n, x)) == std::isnan(x));
        }
    }
}

TEST_CASE("eval_hermite preserves high-order accuracy", "[eval_hermite][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_orthogonal_eval.py#L238-L241
    constexpr double expected = -1.457076485701412e60;
    const double out = xsf::eval_hermite(70, 1.0);
    REQUIRE(xsf::extended_absolute_error(out, expected) <= 1e-14 * std::abs(expected));
}

TEST_CASE("eval_genlaguerre matches constructed polynomials", "[eval_genlaguerre][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_orthogonal_eval.py#L126-L128
    check_poly(
        [](int n, const std::vector<double> &params) { return laguerre_coefficients(n, params[0]); },
        [](int n, const std::vector<double> &params, double x) { return xsf::eval_genlaguerre(n, params[0], x); },
        {{-0.99, 10.0}}, {0.0, 100.0}, 1e-8
    );
}

TEST_CASE("eval_genlaguerre recurrence overload", "[eval_genlaguerre][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_orthogonal_eval.py#L230-L232
    check_recurrence(
        [](int n, const std::vector<double> &params, double x) { return xsf::eval_genlaguerre(n, params[0], x); },
        [](double n, const std::vector<double> &params, double x) { return xsf::eval_genlaguerre(n, params[0], x); },
        {{-0.99, 10.0}}, {0.0, 100.0}
    );
}

TEST_CASE("eval_genlaguerre supports complex inputs", "[eval_genlaguerre][xsf_tests]") {
    const std::vector<std::complex<double>> xs = {{-0.75, -0.25}, {0.0, -0.5}, {0.5, 0.25}, {2.0, -1.0}};
    check_complex_poly(
        [](int n, const std::vector<double> &params) { return laguerre_coefficients(n, params[0]); },
        [](double n, const std::vector<double> &params, std::complex<double> x) {
            return xsf::eval_genlaguerre(n, params[0], x);
        },
        {{-0.99, 10.0}}, xs, 1e-8
    );
}

TEST_CASE("Laguerre evaluators handle domain and NaN inputs", "[eval_genlaguerre][eval_laguerre][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_orthogonal_eval.py#L23-L26
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_orthogonal_eval.py#L258-L265
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (double alpha : {-2.0, -1.0}) {
        REQUIRE(std::isnan(xsf::eval_genlaguerre(0, alpha, 0.0)));
        REQUIRE(std::isnan(xsf::eval_genlaguerre(0.1, alpha, 0.0)));
        REQUIRE(std::isnan(xsf::eval_genlaguerre(0.1, alpha, std::complex<double>{0.0, 1.0}).real()));
    }
    for (double alpha : {1.0, nan}) {
        for (double x : {2.0, nan}) {
            const bool expected = std::isnan(alpha) || std::isnan(x);
            for (int n : {-1, 0, 1, 2}) {
                CAPTURE(n, alpha, x);
                REQUIRE(std::isnan(xsf::eval_genlaguerre(n, alpha, x)) == expected);
            }
            for (double n : {0.0, 1.0, 2.0, 3.2}) {
                CAPTURE(n, alpha, x);
                REQUIRE(std::isnan(xsf::eval_genlaguerre(n, alpha, x)) == expected);
                REQUIRE(std::isnan(xsf::eval_genlaguerre(n, alpha, std::complex<double>{x, 0.0}).real()) == expected);
            }
        }
    }
    for (int n : {-1, 0, 1, 2}) {
        REQUIRE(std::isnan(xsf::eval_laguerre(n, nan)));
        REQUIRE(std::isnan(xsf::eval_laguerre(static_cast<double>(n), nan)));
    }
    REQUIRE(xsf::eval_genlaguerre(-1, 0.5, 2.0) == 0.0);
    REQUIRE(xsf::eval_laguerre(-1, 2.0) == 0.0);
    REQUIRE(xsf::eval_laguerre(1, 1.0) == 0.0);
}

TEST_CASE("Laguerre evaluators support noninteger degrees", "[eval_genlaguerre][eval_laguerre][xsf_tests]") {
    // # Generate references with mpmath
    // import random
    // import mpmath as mp
    // mp.mp.dps = 50
    // rng = random.Random(42)
    // for _ in range(10):
    //     n, alpha = rng.uniform(-0.5, 5), rng.uniform(-0.5, 2)
    //     x = rng.uniform(-2, 2)
    //     z = complex(rng.uniform(-2, 2), rng.uniform(-2, 2))
    //     real = [mp.laguerre(n, a, x) for a in (alpha, 0)]
    //     comp = [mp.laguerre(n, a, z) for a in (alpha, 0)]
    //     print(n, alpha, x, z, *map(float, real), *map(complex, comp))
    using test_case = std::tuple<
        double, double, double, std::complex<double>, double, double, std::complex<double>, std::complex<double>>;
    auto [n, alpha, x, z, gen_expected, expected, gen_complex_expected, complex_expected] = GENERATE(
        test_case{
            3.0168473915183607,
            -0.43747311194333266,
            -0.899882726523523,
            {-1.107157047404709, 0.9458848566560496},
            3.3571359978327036,
            5.071956274504558,
            {2.7501488912278984, -5.067058216930386},
            {4.565552958906825, -6.492201083052953}
        },
        test_case{
            3.2218471808260123,
            1.7304489192621135,
            -1.6522446694823354,
            {-0.3123127212589183, -1.8808111222477186},
            35.43987061763628,
            12.328541041551668,
            {0.5545590622098308, 21.40878475834554},
            {-4.914577906462691, 6.650558536076848}
        },
        test_case{
            0.7025088614198185,
            0.7633882202584059,
            -1.8938561212645455,
            {-1.204649397253406, 0.5995377511180928},
            2.618801603517271,
            2.1831307232051684,
            {2.2574036792531826, -0.3297381694009886},
            {1.7932447198143076, -0.36076486040106087}
        },
        test_case{
            2.497178143317692,
            0.05110155510174175,
            0.3570627355036349,
            {1.2377218267113066, -1.974004961287756},
            0.2829141374228878,
            0.22513917244388124,
            {-3.6714650984323414, 0.520271143964418},
            {-3.6256454163053053, 0.38946594339500257}
        },
        test_case{
            3.932005885080444,
            1.2453484874705671,
            -0.6389979339280325,
            {-1.3780820007528738, 1.8288522888271248},
            16.451903856048844,
            4.85681308204293,
            {5.321223956486376, -50.730565429566056},
            {-5.5235047845166605, -23.861650528848717}
        },
        test_case{
            1.3512699981194471,
            -0.26813539154963023,
            -1.613134492666144,
            {1.3899774653898391, 0.4149041254675643},
            3.0428418442147622,
            3.457758705467502,
            {-0.830690356929397, -0.3499319216881307},
            {-0.6509588097953647, -0.3995457877621184}
        },
        test_case{
            3.9392055030090907,
            1.3243294667345447,
            0.14491236581880296,
            {1.8924630559174824, -0.4858624911665861},
            5.617240003306438,
            0.4880626786289751,
            {-1.9789742944141506, -0.31816060808060764},
            {0.12390748939781761, -0.6778114863282778}
        },
        test_case{
            2.5362234720027486,
            1.5735116606324873,
            0.4740790094569842,
            {1.446827601243109, 0.30940858102704816},
            3.3880086090370556,
            0.010277695053089116,
            {-0.14572271889155777, -0.7704139553886191},
            {-0.8811183767068022, -0.030384884119784048}
        },
        test_case{
            3.3751450991820793,
            -0.38543904086084446,
            -1.0884068973938126,
            {-0.8424481455915713, -1.68083209230549},
            5.2701572505753616,
            7.452302216322751,
            {-3.1638821884508115, 8.145531279206233},
            {-2.427851598334995, 10.907452711804467}
        },
        test_case{
            0.7803498749856659,
            -0.2474964264756772,
            -0.8881055875596315,
            {0.5427377770576007, -0.540671284119663},
            1.470851935394412,
            1.6628535580860702,
            {0.34550447723183286, 0.4719608893576973},
            {0.5785111749654909, 0.4488984442803201}
        }
    );
    CAPTURE(n, alpha, x, z);
    REQUIRE(xsf::extended_relative_error(xsf::eval_genlaguerre(n, alpha, x), gen_expected) < 1e-12);
    REQUIRE(xsf::extended_relative_error(xsf::eval_laguerre(n, x), expected) < 1e-12);
    REQUIRE(xsf::extended_relative_error(xsf::eval_genlaguerre(n, alpha, z), gen_complex_expected) < 1e-12);
    REQUIRE(xsf::extended_relative_error(xsf::eval_laguerre(n, z), complex_expected) < 1e-12);
}
