#include "../testing_utils.h"
#include <xsf/tools.h>

struct parabola_root_functor {
    XSF_HOST_DEVICE std::pair<double, double> operator()(double x) const { return {x * x - 1, 2 * x}; }
};

TEST_CASE("newton_root_finder happy path", "[newton_root_finder][xsf_tests]") {
    auto [result, status] = xsf::detail::find_root_newton(parabola_root_functor{}, 2.0);
    CAPTURE(result, status);
    REQUIRE(result == 1.0);
    REQUIRE(status == xsf::detail::NewtonRootFinderStatus::SUCCESS);
}

TEST_CASE("newton_root_finder initial guess is root", "[newton_root_finder][xsf_tests]") {
    auto [result, status] = xsf::detail::find_root_newton(parabola_root_functor{}, 1.0);
    CAPTURE(result, status);
    REQUIRE(result == 1.0);
    REQUIRE(status == xsf::detail::NewtonRootFinderStatus::SUCCESS);
}

TEST_CASE("newton_root_finder respects rtol parameter", "[newton_root_finder][xsf_tests]") {
    auto [result_default_rtol, status_default_rtol] = xsf::detail::find_root_newton(parabola_root_functor{}, 2.0);
    auto [result_loose_rtol, status_loose_rtol] = xsf::detail::find_root_newton(parabola_root_functor{}, 2.0, 1e-4);
    REQUIRE(status_loose_rtol == xsf::detail::NewtonRootFinderStatus::SUCCESS);
    const double rel_error_default_rtol = xsf::extended_relative_error(result_default_rtol, 1.0);
    const double rel_error_loose_rtol = xsf::extended_relative_error(result_loose_rtol, 1.0);
    CAPTURE(
        rel_error_default_rtol, rel_error_loose_rtol, result_default_rtol, status_default_rtol, result_loose_rtol,
        status_loose_rtol
    );
    REQUIRE(rel_error_default_rtol < rel_error_loose_rtol);
}

TEST_CASE("newton_root_finder respects atol parameter", "[newton_root_finder][xsf_tests]") {
    double atol = 1e-4;
    auto [result_default_atol, status_default_atol] = xsf::detail::find_root_newton(parabola_root_functor{}, 2.0);
    auto [result_loose_atol, status_loose_atol] =
        xsf::detail::find_root_newton(parabola_root_functor{}, 2.0, 4 * std::numeric_limits<double>::epsilon(), atol);
    REQUIRE(result_loose_atol < 1.0 + atol);
    REQUIRE(status_loose_atol == xsf::detail::NewtonRootFinderStatus::SUCCESS);
    const double rel_error_default_atol = xsf::extended_relative_error(result_default_atol, 1.0);
    const double rel_error_loose_atol = xsf::extended_relative_error(result_loose_atol, 1.0);
    CAPTURE(
        rel_error_default_atol, rel_error_loose_atol, result_default_atol, status_default_atol, result_loose_atol,
        status_loose_atol
    );
    REQUIRE(rel_error_default_atol < rel_error_loose_atol);
}

TEST_CASE("newton_root_finder max iterations exceeded", "[newton_root_finder][xsf_tests]") {
    auto [result, status] = xsf::detail::find_root_newton(parabola_root_functor{}, 2.0, 1e-14, 0.0, 5);
    CAPTURE(result, status);
    REQUIRE(status == xsf::detail::NewtonRootFinderStatus::MAX_ITERATIONS_EXCEEDED);
}

TEST_CASE("newton_root_finder derivative zero", "[newton_root_finder][xsf_tests]") {
    // At x = 0, the derivative is zero, so the root finder should return NaN and status 2.
    auto [result, status] = xsf::detail::find_root_newton(parabola_root_functor{}, 0.0);
    CAPTURE(result, status);
    REQUIRE(status == xsf::detail::NewtonRootFinderStatus::DERIVATIVE_ZERO);
    REQUIRE(std::isnan(result));
}

TEST_CASE("newton_root_finder objective returns NaN", "[newton_root_finder][xsf_tests]") {
    struct nan_root_functor {
        XSF_HOST_DEVICE std::pair<double, double> operator()(double x) const {
            return {std::numeric_limits<double>::quiet_NaN(), 2 * x};
        }
    };
    auto [result, status] = xsf::detail::find_root_newton(nan_root_functor{}, 2.0);
    CAPTURE(result, status);
    REQUIRE(status == xsf::detail::NewtonRootFinderStatus::OBJECTIVE_RETURNED_NAN);
    REQUIRE(std::isnan(result));
}

TEST_CASE("newton_root_finder objective returns infinity", "[newton_root_finder][xsf_tests]") {
    const double inf = GENERATE(std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());

    struct inf_root_functor {
        double value;

        XSF_HOST_DEVICE std::pair<double, double> operator()(double x) const { return {value, 2 * x}; }
    };

    auto [result, status] = xsf::detail::find_root_newton(inf_root_functor{inf}, 2.0);

    CAPTURE(inf, result, status);
    REQUIRE(status == xsf::detail::NewtonRootFinderStatus::OBJECTIVE_RETURNED_INF);
    REQUIRE(std::isnan(result));
}
