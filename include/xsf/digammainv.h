#pragma once

#include "xsf/cephes/const.h"
#include "xsf/cephes/psi.h"
#include "xsf/cephes/zeta.h"
#include "xsf/config.h"
#include "xsf/error.h"

namespace xsf {
namespace detail {

    // Initial guess for digammainv using Minka (2000):
    //   x0 = exp(y) + 0.5      if y >= -2.22
    //   x0 = -1 / (y + gamma)  if y < -2.22
    // where gamma = 0.5772... is the Euler-Mascheroni constant.
    XSF_HOST_DEVICE inline double digammainv_initial_guess(double y) {
        constexpr double euler_mascheroni = 0.5772156649015328606;
        if (y >= -2.22) {
            return std::exp(y) + 0.5;
        }
        return -1.0 / (y + euler_mascheroni);
    }

    // Refine an initial guess x using Newton-Raphson iteration:
    //   x_new = x_old - (psi(x_old) - y) / psi'(x_old)
    // where psi'(x) = zeta(2, x) (the trigamma function).
    XSF_HOST_DEVICE inline double digammainv_newton_raphson(double x, double y) {
        const int max_iterations = 20; // avoid infinite loops in pathological cases
        for (int i = 0; i < max_iterations; ++i) {
            double psi_x = cephes::psi(x);
            double trigamma_x = cephes::zeta(2, x);
            double step = (psi_x - y) / trigamma_x;
            x -= step;
            if (std::abs(step) <= std::numeric_limits<double>::epsilon() * std::abs(x)) {
                break;
            }
        }
        return x;
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
            return std::numeric_limits<double>::quiet_NaN();
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
