#include "../testing_utils.h"

#include <xsf/wright_bessel.h>

/* Regression tests for gh-291 and gh-292: the factored linear path in
 * wright_bessel_integral overflows once the contour radius asked for by the
 * fit is restored, and the 50-node quadrature can lose all accuracy there.
 * Reference values were computed with mpmath at 150 digits from the defining
 * series Phi(a,b,x) = sum_k x^k / (k! * Gamma(a k + b)) (all terms positive
 * for x > 0, so direct log-space summation is stable). */
TEST_CASE("wright_bessel logspace quadrature gh-291 gh-292", "[wright_bessel][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double>;
    auto [a, b, x, desired_log] = GENERATE(
        // gh-291 headline case: the capped linear path returned 2.469e+189
        // here, 100 orders of magnitude off.
        test_case{1.5, 50.0, 990000.0, 205.99994406954433},
        // gh-292 headline case: returned -inf / NaN before the fix.
        test_case{0.1, 50.0, 3000.0, 1701.3771505224278},
        // exp_term > log(DBL_MAX) takes the log-space path even below the
        // radius cap of 150.
        test_case{0.1, 50.0, 2000.0, 1111.5226928566246},
        // eps_fit > 150 takes the log-space path.
        test_case{2.0, 50.0, 1e7, 126.83444300823496}
    );

    const double log_out = xsf::log_wright_bessel(a, b, x);
    const double log_err = xsf::extended_relative_error(log_out, desired_log);
    CAPTURE(a, b, x, log_out, desired_log, log_err);
    REQUIRE(log_err <= 1e-8);

    if (desired_log < xsf::cephes::detail::MAXLOG) {
        const double out = xsf::wright_bessel(a, b, x);
        const double err = xsf::extended_relative_error(out, std::exp(desired_log));
        CAPTURE(out, err);
        REQUIRE(err <= 1e-6);
    } else {
        // The true value overflows a double: the result must be +inf, never
        // -inf (gh-292).
        const double out = xsf::wright_bessel(a, b, x);
        CAPTURE(out);
        REQUIRE(out == std::numeric_limits<double>::infinity());
    }
}

TEST_CASE("wright_bessel quadrature failure is honest gh-292", "[wright_bessel][xsf_tests]") {
    /* At this point the log-space quadrature has no positive total and even
     * the largest single series term (a rigorous lower bound for the positive
     * Wright function) exceeds DBL_MAX, so the true value provably overflows:
     * +inf for wright_bessel. The logarithm itself cannot be certified by the
     * quadrature here and is reported as NaN instead of a wrong value. */
    const double out = xsf::wright_bessel(0.05, 40.0, 3700.0);
    const double log_out = xsf::log_wright_bessel(0.05, 40.0, 3700.0);
    CAPTURE(out, log_out);
    REQUIRE(out == std::numeric_limits<double>::infinity());
    REQUIRE(std::isnan(log_out));
}
