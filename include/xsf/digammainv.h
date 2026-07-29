#pragma once

#include "xsf/cephes/const.h"
#include "xsf/cephes/psi.h"
#include "xsf/cephes/zeta.h"
#include "xsf/config.h"
#include "xsf/error.h"
#include "xsf/tools.h"

namespace xsf {
namespace detail {

    struct digammainv_root_functor {
        double y;
        XSF_HOST_DEVICE std::pair<double, double> operator()(double x) const {
            return {cephes::psi(x) - y, cephes::zeta(2, x)};
        }
    };

    // Initial guess for digammainv using Minka (2000):
    //   x0 = exp(y) + 0.5      if y >= -2.22
    //   x0 = -1 / (y + gamma)  if y < -2.22
    // where gamma = 0.5772... is the Euler-Mascheroni constant.
    XSF_HOST_DEVICE inline double digammainv_initial_guess(double y) {
        if (y >= -2.22) {
            return std::exp(y) + 0.5;
        }
        return -1.0 / (y + cephes::detail::SCIPY_EULER);
    }

    // Refine an initial guess x using Newton-Raphson iteration:
    //   x_new = x_old - (psi(x_old) - y) / psi'(x_old)
    // where psi'(x) = zeta(2, x) (the trigamma function).
    XSF_HOST_DEVICE inline double digammainv_newton_raphson(double x, double y) {
        auto [res, status] = find_root_newton(digammainv_root_functor{y}, x);
        if (status == NewtonRootFinderStatus::MAX_ITERATIONS_EXCEEDED) {
            set_error("digammainv", SF_ERROR_NO_RESULT, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (status == NewtonRootFinderStatus::DERIVATIVE_ZERO) {
            set_error("digammainv", SF_ERROR_SINGULAR, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (status == NewtonRootFinderStatus::OBJECTIVE_RETURNED_NAN) {
            set_error("digammainv", SF_ERROR_NO_RESULT, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (status == NewtonRootFinderStatus::OBJECTIVE_RETURNED_INF) {
            set_error("digammainv", SF_ERROR_OVERFLOW, NULL);
            return res;
        }
        return res;
    }

    // Compute the inverse of the digamma function for real double input y,
    // i.e. find x > 0 such that psi(x) = y.
    XSF_HOST_DEVICE inline double digammainv(double y) {
        // Reference:
        //   T. Minka, "Estimating a Dirichlet distribution", 2000. Appendix C
        //   https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
        if (std::isnan(y)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        // digamma(x) app. log(x) for large x, so digammainv(y) must be inf for y > log(DBL_MAX)
        if (y > cephes::detail::MAXLOG) {
            set_error("digammainv", SF_ERROR_OVERFLOW, NULL);
            return std::numeric_limits<double>::infinity();
        }
        // digamma(x) -> - inf for x -> 0+, so digammainv must be zero for y = - inf
        if (std::isinf(y) && y < 0) {
            return 0.0;
        }
        double x = digammainv_initial_guess(y);
        return digammainv_newton_raphson(x, y);
    }

} // namespace detail

XSF_HOST_DEVICE inline double digammainv(double y) { return detail::digammainv(y); }

XSF_HOST_DEVICE inline float digammainv(float y) {
    return static_cast<float>(detail::digammainv(static_cast<double>(y)));
}

} // namespace xsf
