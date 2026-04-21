#pragma once

#include <cmath>
#include <limits>

#include "xsf/cephes/psi.h"
#include "xsf/cephes/zeta.h"

namespace xsf {
namespace detail {

    // Initial guess for digamma_inv using Minka (2000):
    //   x0 = exp(y) + 0.5      if y >= -2.22
    //   x0 = -1 / (y + gamma)  if y < -2.22
    // where gamma = 0.5772... is the Euler-Mascheroni constant.
    inline double digamma_inv_initial_guess(double y) {
        constexpr double euler_mascheroni = 0.5772156649015328606;
        if (y >= -2.22) {
            return std::exp(y) + 0.5;
        }
        return -1.0 / (y + euler_mascheroni);
    }

    // Refine an initial guess x using Newton-Raphson iteration:
    //   x_new = x_old - (psi(x_old) - y) / psi'(x_old)
    // where psi'(x) = zeta(2, x) (the trigamma function).
    // x is passed by value; the refined result is returned.
    inline double digamma_inv_newton_raphson(double x, double y) {
        const int max_iterations = 20; // avoid infinite loops in pathological cases
        for (int i = 0; i < max_iterations; ++i) {
            double psi_x = cephes::psi(x);
            double psi1_x = cephes::zeta(2, x); // trigamma: psi'(x) = zeta(2, x)
            double step = (psi_x - y) / psi1_x;
            x -= step;
            if (std::abs(step) <= std::numeric_limits<double>::epsilon() * std::abs(x)) {
                break;
            }
        }
        return x;
    }

    // Compute the inverse of the digamma function for real double input y,
    // i.e. find x > 0 such that psi(x) = y.
    inline double digamma_inv(double y) {
        // Reference:
        //   T. Minka, "Estimating a Dirichlet distribution", 2000. Appendix C
        //   https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
        if (std::isnan(y)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (std::isinf(y)) {
            // digamma(x) -> +inf as x -> +inf, and digamma(x) -> -inf as x -> 0+.
            return y > 0 ? std::numeric_limits<double>::infinity() : 0.0;
        }

        double x = digamma_inv_initial_guess(y);
        return digamma_inv_newton_raphson(x, y);
    }

} // namespace detail

inline double digamma_inv(double y) { return detail::digamma_inv(y); }

inline float digamma_inv(float y) { return static_cast<float>(detail::digamma_inv(static_cast<double>(y))); }

} // namespace xsf
