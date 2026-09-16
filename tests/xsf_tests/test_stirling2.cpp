#include "../testing_utils.h"
#include <xsf/stirling2.h>

// Known values of S(n, k) for n = 0..10 and k = 0..n.
// Row n contains S(n, 0), S(n, 1), ..., S(n, n).
static const double stirling2_table[11][11] = {
    {1},
    {0, 1},
    {0, 1, 1},
    {0, 1, 3, 1},
    {0, 1, 7, 6, 1},
    {0, 1, 15, 25, 10, 1},
    {0, 1, 31, 90, 65, 15, 1},
    {0, 1, 63, 301, 350, 140, 21, 1},
    {0, 1, 127, 966, 1701, 1050, 266, 28, 1},
    {0, 1, 255, 3025, 7770, 6951, 2646, 462, 36, 1},
    {0, 1, 511, 9330, 34105, 42525, 22827, 5880, 750, 45, 1},
};

TEST_CASE("stirling2 table values", "[stirling2][xsf_tests]") {
    const double rtol = 1e-12;
    const float rtol_f = 1e-6f;
    for (int n = 0; n <= 10; ++n) {
        for (int k = 0; k <= n; ++k) {
            const double result = xsf::stirling2(static_cast<double>(n), static_cast<double>(k));
            const double expected = stirling2_table[n][k];
            const auto rel_error = xsf::extended_relative_error(result, expected);
            CAPTURE(n, k, result, expected, rel_error);
            REQUIRE(rel_error <= rtol);

            const float n_f = static_cast<float>(n);
            const float k_f = static_cast<float>(k);
            const float result_f = xsf::stirling2(n_f, k_f);
            const float expected_f = static_cast<float>(expected);
            const auto rel_error_f = xsf::extended_relative_error(result_f, expected_f);
            CAPTURE(n_f, k_f, result_f, expected_f, rel_error_f);
            REQUIRE(rel_error_f <= rtol_f);
        }
    }
}

TEST_CASE("stirling2 single values", "[stirling2][xsf_tests]") {
    const double rtol = 1e-12;

    auto check = [&](double n, double k, double expected) {
        const double result = xsf::stirling2(n, k);
        const auto rel_error = xsf::extended_relative_error(result, expected);
        CAPTURE(n, k, result, expected, rel_error);
        REQUIRE(rel_error <= rtol);
    };

    check(0.0, 0.0, stirling2_table[0][0]);
    check(4.0, 2.0, stirling2_table[4][2]);
    check(5.0, 3.0, 25.0);
}

TEST_CASE("stirling2 negative inputs", "[stirling2][xsf_tests]") {
    REQUIRE(xsf::stirling2(-1.0, -1.0) == 0.0);
    REQUIRE(xsf::stirling2(-1.0, 2.0) == 0.0);
    REQUIRE(xsf::stirling2(2.0, -1.0) == 0.0);
}

TEST_CASE("stirling2 mixed inputs", "[stirling2][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
    // n=-1 or k=-2 returns 0; otherwise matches known table entry.
    auto [n, k, expected] = GENERATE(
        test_case{-1.0, -2.0, 0.0}, test_case{0.0, 0.0, 1.0}, test_case{3.0, 2.0, 3.0}, test_case{5.0, 3.0, 25.0},
        test_case{8.0, 5.0, 1050.0}, test_case{10.0, 7.0, 5880.0}, test_case{10.0, 3.0, 9330.0}
    );
    const double rtol = 1e-13;
    const double result = xsf::stirling2(n, k);
    const auto rel_error = xsf::extended_relative_error(result, expected);
    CAPTURE(n, k, result, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("stirling2 k out of range", "[stirling2][xsf_tests]") {
    // k > n or k = 0 with n > 0 gives 0
    REQUIRE(xsf::stirling2(3.0, 4.0) == 0.0);
    REQUIRE(xsf::stirling2(5.0, 0.0) == 0.0);
    REQUIRE(xsf::stirling2(0.0, 1.0) == 0.0);
}

TEST_CASE("stirling2 temme approximation accuracy", "[stirling2][xsf_tests]") {
    // For n > 50 the Temme approximation is used. We verify it against
    // analytically exact identities that hold for all n:
    //   S(n, 1)   = 1                  (exactly one way to put all in one subset)
    //   S(n, n)   = 1                  (each element is its own subset)
    //   S(n, n-1) = n*(n-1)/2 = C(n,2) (choose which two elements share a subset)
    //   S(n, 2)   = 2^(n-1) - 1        (exact for n <= 53 in double)
    const double max_rtol = 2e-5;

    for (double n : {51.0, 55.0, 60.0, 75.0, 100.0}) {
        // S(n, 1) = 1
        {
            const double result = xsf::stirling2(n, 1.0);
            const auto rel_error = xsf::extended_relative_error(result, 1.0);
            CAPTURE(n, 1.0, result, rel_error);
            REQUIRE(rel_error <= max_rtol);
        }

        // S(n, n) = 1
        {
            const double result = xsf::stirling2(n, n);
            const auto rel_error = xsf::extended_relative_error(result, 1.0);
            CAPTURE(n, n, result, rel_error);
            REQUIRE(rel_error <= max_rtol);
        }

        // S(n, n-1) = n*(n-1)/2
        {
            const double expected = n * (n - 1.0) / 2.0;
            const double result = xsf::stirling2(n, n - 1.0);
            const auto rel_error = xsf::extended_relative_error(result, expected);
            CAPTURE(n, n - 1.0, result, expected, rel_error);
            REQUIRE(rel_error <= max_rtol);
        }

        // S(n, 2) = 2^(n-1) - 1  (only test where representable in double, n <= 53)
        if (n <= 53.0) {
            const double expected = std::pow(2.0, n - 1.0) - 1.0;
            const double result = xsf::stirling2(n, 2.0);
            const auto rel_error = xsf::extended_relative_error(result, expected);
            CAPTURE(n, 2.0, result, expected, rel_error);
            REQUIRE(rel_error <= max_rtol);
        }
    }
}
