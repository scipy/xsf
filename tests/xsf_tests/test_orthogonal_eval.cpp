#include "../testing_utils.h"

#include <xsf/orthogonal_eval.h>

TEST_CASE("eval_jacobi real values", "[eval_jacobi][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double, double>;

    auto [n, alpha, beta, x, expected] = GENERATE(
        test_case{2, 0.0, 0.0, 0.5, -0.125}, test_case{3, 1.0, 1.0, -0.3, 0.71100000000000019},
        test_case{4, 0.5, 1.5, 0.2, 0.54993750000000008}, test_case{0, 2.0, 3.0, 0.7, 1.0},
        test_case{1, -0.5, 0.5, 0.4, -0.099999999999999978}, test_case{5, 2.0, -0.5, -0.6, -0.27551999999999993}
    );

    const double result = xsf::eval_jacobi(n, alpha, beta, x);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(n, alpha, beta, x, result, expected, error);
    REQUIRE(error <= 1e-12);
}

TEST_CASE("eval_jacobi_l integer order", "[eval_jacobi][xsf_tests]") {
    using test_case = std::tuple<int, double, double, double, double>;

    auto [n, alpha, beta, x, expected] = GENERATE(
        test_case{2, 0.0, 0.0, 0.5, -0.125}, test_case{3, 1.0, 1.0, -0.3, 0.71100000000000019},
        test_case{4, 0.5, 1.5, 0.2, 0.54993750000000008}, test_case{0, 2.0, 3.0, 0.7, 1.0},
        test_case{1, -0.5, 0.5, 0.4, -0.099999999999999978}, test_case{5, 2.0, -0.5, -0.6, -0.27551999999999988}
    );

    const double result = xsf::eval_jacobi_l(n, alpha, beta, x);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(n, alpha, beta, x, result, expected, error);
    REQUIRE(error <= 1e-12);
}

TEST_CASE("eval_jacobi degenerate alpha=-1 |beta|=1", "[eval_jacobi][xsf_tests]") {
    using test_case = std::tuple<int, double, double>;

    auto [n, x, expected] = GENERATE(
        test_case{0, 0.3, 1.0}, test_case{1, 0.3, -0.69999999999999996}, test_case{2, 0.3, -0.315},
        test_case{3, 0.3, 0.19250000000000023}
    );

    const double from_double = xsf::eval_jacobi(static_cast<double>(n), -1.0, 1.0, x);
    const double from_l = xsf::eval_jacobi_l(n, -1.0, 1.0, x);
    CAPTURE(n, x, from_double, from_l, expected);
    REQUIRE(xsf::extended_relative_error(from_double, expected) <= 1e-12);
    REQUIRE(xsf::extended_relative_error(from_l, expected) <= 1e-12);
}

TEST_CASE("eval_jacobi float overloads", "[eval_jacobi][xsf_tests]") {
    auto [n, alpha, beta, x] = GENERATE(
        std::tuple<float, float, float, float>{2.0F, 0.0F, 0.0F, 0.5F},
        std::tuple<float, float, float, float>{4.0F, 0.5F, 1.5F, 0.2F},
        std::tuple<float, float, float, float>{5.0F, 2.0F, -0.5F, -0.6F}
    );

    const float result = xsf::eval_jacobi(n, alpha, beta, x);
    const float expected = static_cast<float>(xsf::eval_jacobi(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(beta), static_cast<double>(x)
    ));
    CAPTURE(n, alpha, beta, x, result, expected);
    REQUIRE(result == expected);
}

TEST_CASE("eval_sh_jacobi real values", "[eval_sh_jacobi][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double, double>;

    auto [n, p, q, x, expected] = GENERATE(
        test_case{2, 3.0, 1.0, 0.5, -0.016666666666666698}, test_case{3, 2.0, 1.0, 0.3, 0.01128571428571429},
        test_case{4, 4.0, 2.0, 0.7, -0.0035363636363636357}, test_case{0, 2.0, 1.0, 0.4, 1.0},
        test_case{1, 3.0, 2.0, 0.6, 0.099999999999999978}, test_case{5, 5.0, 3.0, 0.2, 0.00079552447552447498}
    );

    const double result = xsf::eval_sh_jacobi(n, p, q, x);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(n, p, q, x, result, expected, error);
    REQUIRE(error <= 1e-12);
}

TEST_CASE("eval_sh_jacobi_l integer order", "[eval_sh_jacobi][xsf_tests]") {
    using test_case = std::tuple<int, double, double, double, double>;

    auto [n, p, q, x, expected] = GENERATE(
        test_case{2, 3.0, 1.0, 0.5, -0.016666666666666673}, test_case{3, 2.0, 1.0, 0.3, 0.011285714285714286},
        test_case{4, 4.0, 2.0, 0.7, -0.0035363636363636352}, test_case{0, 2.0, 1.0, 0.4, 1.0},
        test_case{1, 3.0, 2.0, 0.6, 0.099999999999999978}, test_case{5, 5.0, 3.0, 0.2, 0.00079552447552447584}
    );

    const double result = xsf::eval_sh_jacobi_l(n, p, q, x);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(n, p, q, x, result, expected, error);
    REQUIRE(error <= 1e-12);
}

TEST_CASE("eval_sh_jacobi float overloads", "[eval_sh_jacobi][xsf_tests]") {
    auto [n, p, q, x] = GENERATE(
        std::tuple<float, float, float, float>{2.0F, 3.0F, 1.0F, 0.5F},
        std::tuple<float, float, float, float>{4.0F, 4.0F, 2.0F, 0.7F}
    );

    const float result = xsf::eval_sh_jacobi(n, p, q, x);
    const float expected = static_cast<float>(xsf::eval_sh_jacobi(
        static_cast<double>(n), static_cast<double>(p), static_cast<double>(q), static_cast<double>(x)
    ));
    CAPTURE(n, p, q, x, result, expected);
    REQUIRE(result == expected);
}
