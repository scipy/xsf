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

namespace xsf {
namespace cephes {

    /* Compute acos(x) / pi */
    template<typename T>
    XSF_HOST_DEVICE T acospi(double x)
    {  
        double r = std::acos(x)/M_PI;
        if (std::isgreater(r, 1.0)) {  
            return 1.0;
        }
        return r;
    }

    /* Compute asin(x) / pi */
    template<typename T>
    XSF_HOST_DEVICE T asinpi(double x)
    {  
        double r = std::asin(x)/M_PI;
        if (std::isgreater(r, 1.0)) {  
            return 1.0;
        }
        return r;
    }

    /* Compute atan(x) / pi */
    template<typename T>
    XSF_HOST_DEVICE T atanpi(T x)
    {
        double r = std::atan(x)/M_PI;
        if (isgreater(std::fabs(r), 0.5)) {
            return std::copysign(0.5, r);
        }
        return r;
    }

    /* Compute atan2(y, x) / pi */
    template<typename T>
    XSF_HOST_DEVICE T atan2pi(double y, double x)
    {
        double r = std::atan2(y, x)/M_PI;
        if (std::isgreater(std::fabs(r), 1.0)) {
            return std::copysign(1.0, r);
        }
        return r;
    }

    /* Compute sin(pi * x). */
    template <typename T>
    XSF_HOST_DEVICE T sinpi(T x) {
        T s = 1.0;

        if (std::signbit(x)) {
            x = -x;
            s = -1.0;
        }

        T r = std::fmod(x, 2.0);
        if (r < 0.5) {
            return s * std::sin(M_PI * r);
        } else if (r > 1.5) {
            return s * std::sin(M_PI * (r - 2.0));
        } else {
            return -s * std::sin(M_PI * (r - 1.0));
        }
    }

    /* Compute cos(pi * x) */
    template <typename T>
    XSF_HOST_DEVICE T cospi(T x) {
        if (x < 0.0) {
            x = -x;
        }

        T r = std::fmod(x, 2.0);
        if (r == 0.5) {
            // We don't want to return -0.0
            return 0.0;
        }
        if (r < 1.0) {
            return -std::sin(M_PI * (r - 0.5));
        } else {
            return std::sin(M_PI * (r - 1.5));
        }
    }

    /* Compute tan(pi * x). */
    template <typename T>
    XSF_HOST_DEVICE T tanpi(T x) {
        T y, absy;         
        if (!std::isfinite(x)) {
            return std::tan(x);
        }
        y = x - 2.0 * std::round(0.5 * x);
        absy = std::fabs(y);
        if (absy == 0.0) {
            return std::copysign(0.0, x);
        } 
        if (absy == 1.0) {
            return std::copysign(0.0, -x);
        }
        if (absy == 0.5) {
            errno = ERANGE;
            return 1.0 / std::copysign(0.0, y);
        }
        if (absy > 0.5) {
            y -= std::copysign(1.0, y);
            absy = std::fabs(y);
        }
        if (absy <= 0.25) {
            return std::tan(M_PI * y);
        }
        return std::copysign(1.0 / std::tan(M_PI * (0.5 - absy)), y);
    }

} // namespace cephes
} // namespace xsf
