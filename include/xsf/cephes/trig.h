/* Translated into C++ by SciPy developers in 2024.
 *
 * Original author: Josh Wilson, 2020.
 */

/*
 * Implement half-turn trig functions for real x. Since the periods
 * of these functions are integral (and thus representable in double
 * precision), it's possible to compute e.g., sinpi(x) them with
 * greater accuracy than sin(pi * x).
 */
#pragma once

#include "../config.h"
#include "numbers.h"

namespace xsf {
namespace cephes {

    /* Compute acos(x) / pi */
    template <typename T>
    XSF_HOST_DEVICE T acospi(double x) {
        double r = std::acos(x) / pi_v<T>;
        if (std::isgreater(r, T(1.0))) {
            return T(1.0);
        }
        return r;
    }

    /* Compute asin(x) / pi */
    template <typename T>
    XSF_HOST_DEVICE T asinpi(double x) {
        double r = std::asin(x) / pi_v<T>;
        if (std::isgreater(r, T(1.0))) {
            return T(1.0);
        }
        return r;
    }

    /* Compute atan(x) / pi */
    template <typename T>
    XSF_HOST_DEVICE T atanpi(T x) {
        double r = std::atan(x) / pi_v<T>;
        if (isgreater(std::fabs(r), T(0.5))) {
            return std::copysign(T(0.5), r);
        }
        return r;
    }

    /* Compute atan2(y, x) / pi */
    template <typename T>
    XSF_HOST_DEVICE T atan2pi(double y, double x) {
        double r = std::atan2(y, x) / pi_v<T>;
        if (std::isgreater(std::fabs(r), T(1.0))) {
            return std::copysign(T(1.0), r);
        }
        return r;
    }

    /* Compute sin(pi * x). */
    template <typename T>
    XSF_HOST_DEVICE T sinpi(T x) {
        T s = T(1.0);

        if (std::signbit(x)) {
            x = -x;
            s = -T(1.0);
        }

        T r = std::fmod(x, T(2.0));
        if (r < T(0.5)) {
            return s * std::sin(pi_v<T> * r);
        } else if (r > T(1.5)) {
            return s * std::sin(pi_v<T> * (r - T(2.0)));
        } else {
            return -s * std::sin(pi_v<T> * (r - T(1.0)));
        }
    }

    /* Compute cos(pi * x) */
    template <typename T>
    XSF_HOST_DEVICE T cospi(T x) {
        if (x < T(0.0)) {
            x = -x;
        }

        T r = std::fmod(x, T(2.0));
        if (r == T(0.5)) {
            // We don't want to return -T(0.0)
            return T(0.0);
        }
        if (r < T(1.0)) {
            return -std::sin(pi_v<T> * (r - T(0.5)));
        } else {
            return std::sin(pi_v<T> * (r - T(1.5)));
        }
    }

    /* Compute tan(pi * x). */
    template <typename T>
    XSF_HOST_DEVICE T tanpi(T x) {
        T y, absy;
        if (!std::isfinite(x)) {
            return std::tan(x);
        }
        y = x - T(2.0) * std::round(T(0.5) * x);
        absy = std::fabs(y);
        if (absy == T(0.0)) {
            return std::copysign(T(0.0), x);
        }
        if (absy == T(1.0)) {
            return std::copysign(T(0.0), -x);
        }
        if (absy == T(0.5)) {
            errno = ERANGE;
            return T(1.0) / std::copysign(T(0.0), y);
        }
        if (absy > T(0.5)) {
            y -= std::copysign(T(1.0), y);
            absy = std::fabs(y);
        }
        if (absy <= T(0.25)) {
            return std::tan(pi_v<T> * y);
        }
        return std::copysign(T(1.0) / std::tan(pi_v<T> * (T(0.5) - absy)), y);
    }

} // namespace cephes
} // namespace xsf
