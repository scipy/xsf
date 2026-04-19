#pragma once

#include <cmath>
#include <limits>

#include "cephes/psi.h"
#include "cephes/zeta.h"

namespace xsf {
namespace detail {

    // Compute the inverse of the digamma function for real double input y,
    // i.e. find x > 0 such that psi(x) = y.
    //
    // Algorithm: Newton-Raphson iteration using
    //   x_new = x_old - (psi(x_old) - y) / psi'(x_old)
    // where psi'(x) = zeta(2, x) (the trigamma function).
    //
    // Initial guess (Minka 2000):
    //   x0 = exp(y) + 0.5      if y >= -2.22
    //   x0 = -1 / (y + gamma)  if y < -2.22
    // where gamma = 0.5772... is the Euler-Mascheroni constant.
    //
    // Reference:
    //   T. Minka, "Estimating a Dirichlet distribution", 2000.
    //   https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
    inline double inv_digamma(double y) {
        if (std::isnan(y)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (std::isinf(y)) {
            return y > 0 ? std::numeric_limits<double>::infinity() : 0.0;
        }

        constexpr double euler_mascheroni = 0.5772156649015328606;

        // Initial guess
        double x;
        if (y >= -2.22) {
            x = std::exp(y) + 0.5;
        } else {
            x = -1.0 / (y + euler_mascheroni);
        }

        // Newton-Raphson iterations
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

} // namespace detail

inline double inv_digamma(double y) { return detail::inv_digamma(y); }

inline float inv_digamma(float y) { return static_cast<float>(detail::inv_digamma(static_cast<double>(y))); }

} // namespace xsf
