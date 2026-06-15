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
            // We don't want to return -0.0
            return T(0.0);
        }
        if (r < T(1.0)) {
            return -std::sin(pi_v<T> * (r - T(0.5)));
        } else {
            return std::sin(pi_v<T> * (r - T(1.5)));
        }
    }
} // namespace cephes
} // namespace xsf
