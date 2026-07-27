#include "../../include/xsf/config.h"
#include "../testing_utils.h"

#include <xsf/cephes/polevl.h>
#include <xsf/evalpoly.h>
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

double polyval(const std::vector<double> &coeffs, double x) {
    if (coeffs.size() == 1) {
        return coeffs[0];
    }

    const std::vector<double> reversed(coeffs.rbegin(), coeffs.rend());
    return xsf::cephes::polevl(x, reversed.data(), reversed.size() - 1);
}

std::complex<double> polyval(const std::vector<double> &coeffs, std::complex<double> x) {
    if (coeffs.size() == 1) {
        return coeffs[0];
    }

    const std::vector<double> reversed(coeffs.rbegin(), coeffs.rend());
    return xsf::cevalpoly(reversed.data(), reversed.size() - 1, x);
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
        [](std::ptrdiff_t n, const std::vector<double> &params, double x) {
            return xsf::eval_jacobi(n, params[0], params[1], x);
        },
        [](double n, const std::vector<double> &params, double x) {
            return xsf::eval_jacobi(n, params[0], params[1], x);
        },
        {{-0.99, 10.0}, {-0.99, 10.0}}, {-1.0, 1.0}
    );
}

TEST_CASE("eval_sh_jacobi recurrence overload", "[eval_sh_jacobi][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_orthogonal_eval.py#L190-L192
    check_recurrence(
        [](std::ptrdiff_t n, const std::vector<double> &params, double x) {
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
        auto out = xsf::eval_jacobi(static_cast<std::ptrdiff_t>(n), -1.0, 1.0, x);
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
        auto out = xsf::eval_jacobi(static_cast<std::ptrdiff_t>(n), -1.0, -1.0, x);
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
