#pragma once

#include "config.h"
#include "cephes/ellpk.h"

namespace xsf {

namespace detail {

    XSF_HOST_DEVICE inline double agm_iter(double a, double b) {
        // Arithmetic-geometric mean, iterative implementation
        // a and b must be positive (not zero, not nan).
        double amean = 0.5 * a + 0.5 * b;
        int count = 20;
        while (count > 0 && (amean != a && amean != b)) {
            double gmean = std::sqrt(a) * std::sqrt(b);
            a = amean;
            b = gmean;
            amean = 0.5 * a + 0.5 * b;
            count -= 1;
        }
        return amean;
    }

} // namespace detail

// Arithmetic-geometric mean.
XSF_HOST_DEVICE inline double agm(double a, double b) {
    // sqrthalfmax is std::sqrt(std::numeric_limits<double>::max() / 2)
    // invsqrthalfmax is 1/sqrthalfmax
    constexpr double sqrthalfmax = 9.480751908109176e+153;
    constexpr double invsqrthalfmax = 1.0547686614863e-154;

    if (std::isnan(a) || std::isnan(b)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if ((a < 0 && b > 0) || (a > 0 && b < 0)) {
        // a and b have opposite sign.
        return std::numeric_limits<double>::quiet_NaN();
    }

    if ((std::isinf(a) || std::isinf(b)) && (a == 0 || b == 0)) {
        // One value is inf and the other is 0.
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (a == 0 || b == 0) {
        // At least one of the arguments is 0.
        return 0.0;
    }

    if (a == b) {
        return a;
    }

    int sgn = 1;
    if (a < 0) {
        sgn = -1;
        a = -a;
        b = -b;
    }

    // At this point, a and b are both positive and not nan.

    if ((invsqrthalfmax < a && a < sqrthalfmax) && (invsqrthalfmax < b && b < sqrthalfmax)) {
        double e = 4 * a * b / ((a + b) * (a + b));
        return sgn * (M_PI / 4) * (a + b) / cephes::ellpk(e);
    }

    // At least one value is "extreme" (very big or very small).
    // Use the iteration to avoid overflow or underflow.

    return sgn * detail::agm_iter(a, b);
}

inline float agm(float a, float b) { return agm(static_cast<double>(a), static_cast<double>(b)); }

} // namespace xsf
