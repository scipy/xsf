#include "../testing_utils.h"

#include <xsf/legendre.h>

// backport of std::mdspan (since C++23)
#define MDSPAN_USE_PAREN_OPERATOR 1
#include <xsf/third_party/kokkos/mdspan.hpp>

// From https://github.com/scipy/scipy/blob/bdd3b0e/scipy/special/tests/test_legendre.py#L693-L697
TEST_CASE("lqmn TestLegendreFunctions.test_lqmn", "[lqmn][lqn][real][smoketest]") {
    constexpr double atol = 1.5e-4;

    constexpr int m = 0;
    constexpr int n = 2;
    constexpr double x = 0.5;

    constexpr int m1p = m + 1;
    constexpr int n1p = n + 1;
    // lqmn requires buffer space for at least 2x2
    constexpr int bufsize = std::max(2, m1p) * std::max(2, n1p);

    // lqmnf = special.lqmn(0, 2, .5)
    double lqmnf0_data[bufsize], lqmnf1_data[bufsize];
    auto lqmnf0 = std::mdspan(lqmnf0_data, m1p, n1p);
    auto lqmnf1 = std::mdspan(lqmnf1_data, m1p, n1p);
    xsf::lqmn(x, lqmnf0, lqmnf1);

    // lqf = special.lqn(2, .5)
    double lqf0_data[n1p], lqf1_data[n1p];
    auto lqf0 = std::mdspan(lqf0_data, n1p);
    auto lqf1 = std::mdspan(lqf1_data, n1p);
    xsf::lqn(x, lqf0, lqf1);

    // assert_allclose(lqmnf[0][0], lqf[0], atol=1.5e-4, rtol=0)
    for (int n = 0; n < n1p; ++n) {
        auto error0 = xsf::extended_absolute_error(lqmnf0(0, n), lqf0(n));
        CAPTURE(n, lqmnf0(0, n), lqf0(n), error0);
        REQUIRE(error0 <= atol);
    }

    // assert_allclose(lqmnf[1][0], lqf[1], atol=1.5e-4, rtol=0)
    for (int n = 0; n < n1p; ++n) {
        auto error1 = xsf::extended_absolute_error(lqmnf1(0, n), lqf1(n));
        CAPTURE(n, lqmnf1(0, n), lqf1(n), error1);
        REQUIRE(error1 <= atol);
    }

    // some additional smoke checks
    for (int n = 0; n < n1p; ++n) {
        CAPTURE(n, lqmnf0(0, n));
        REQUIRE(lqmnf0(0, n) != 0.0);
        REQUIRE(std::isfinite(lqmnf0(0, n)));
    }
    for (int n = 0; n < n1p; ++n) {
        CAPTURE(n, lqmnf1(0, n));
        REQUIRE(lqmnf1(0, n) != 0.0);
        REQUIRE(std::isfinite(lqmnf1(0, n)));
    }
}

// From https://github.com/scipy/scipy/blob/bdd3b0e/scipy/special/tests/test_legendre.py#L699-L708
TEST_CASE("lqmn TestLegendreFunctions.test_lqmn_gt1", "[lqmn][real][smoketest]") {
    constexpr double atol = 1.5e-7;

    constexpr int m = 2;
    constexpr int n = 1;
    constexpr int m1p = m + 1;
    constexpr int n1p = n + 1;

    double lqmnf0_data[m1p * n1p], lqmnf1_data[m1p * n1p];
    auto lqmnf0 = std::mdspan(lqmnf0_data, m1p, n1p);
    auto lqmnf1 = std::mdspan(lqmnf1_data, m1p, n1p);

    // algorithm for real arguments changes at 1.0001
    // test against analytical result for m=2, n=1
    auto x0 = 1.0001;
    auto delta = 0.00002;

    // for x in (x0-delta, x0+delta):
    for (double x : {x0 - delta, x0 + delta}) {
        // lq = special.lqmn(2, 1, x)[0][-1, -1]
        xsf::lqmn(x, lqmnf0, lqmnf1);
        double lq = lqmnf0(m, n); // _[-1, -1] corresponds to _[m, n]

        // expected = 2/(x*x-1)
        double expected = 2.0 / (x * x - 1.0);

        // assert_allclose(lq, expected, atol=1.5e-7, rtol=0)
        auto error = xsf::extended_absolute_error(lq, expected);
        CAPTURE(x, lq, expected, error, atol);
        REQUIRE(error <= atol);
    }
}

TEST_CASE("lqmn complex", "[lqmn][complex][smoketest]") {
    constexpr double atol = 1e-16;
    constexpr double x = 0.5;
    // lqmn requires buffer space for at least 2x2
    constexpr int bufsize = 2 * 2;

    // (q_mn, qp_mn) = lqmn(0, 0, 0.5)
    double q_data[bufsize], qp_data[bufsize];
    auto q_mn = std::mdspan(q_data, 1, 1);
    auto qp_mn = std::mdspan(qp_data, 1, 1);
    xsf::lqmn(x, q_mn, qp_mn);
    auto q = q_mn(0, 0);
    auto qp = qp_mn(0, 0);

    // (cq_mn, cqp_mn) = lqmn(0, 0, 0.5 + 0j)
    std::complex<double> cq_data[bufsize], cqp_data[bufsize];
    auto cq_mn = std::mdspan(cq_data, 1, 1);
    auto cqp_mn = std::mdspan(cqp_data, 1, 1);
    xsf::lqmn(std::complex<double>(x, 0.0), cq_mn, cqp_mn);
    auto cq = cq_mn(0, 0);
    auto cqp = cqp_mn(0, 0);

    // abs(_.imag) <= atol
    CAPTURE(q, qp, cq, cqp, atol);
    CHECK(std::abs(std::imag(cq)) <= atol);
    CHECK(std::abs(std::imag(cqp)) <= atol);

    // abs(q - cq.real) <= atol
    auto err_q = xsf::extended_absolute_error(q, std::real(cq));
    CAPTURE(err_q);
    REQUIRE(err_q <= atol);

    // abs(qp - cqp.real) <= atol
    auto err_qp = xsf::extended_absolute_error(qp, std::real(cqp));
    CAPTURE(err_qp);
    REQUIRE(err_qp <= atol);
}
