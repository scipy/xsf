// Translated from Cython to C++ by the xsf developers in 2026.
#pragma once

#include "../cephes/poch.h"
#include "../config.h"
#include "../error.h"
#include "specfun.h"

namespace xsf {

inline double hyperu(double a, double b, double x) {
    if (std::isnan(a) || std::isnan(b) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (x < 0.0) {
        set_error("hyperu", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (x == 0.0) {
        if (b > 1.0) {
            // odd floor(a) gives (-1)^m for a = -m and the sign of 1 / Gamma(a) otherwise
            const double fa = std::floor(a);
            const bool finite = std::isfinite(a);
            const bool odd = a < 0.0 && finite && std::fmod(fa, 2.0) != 0.0;
            if (a <= 0.0 && a == fa && finite) {
                // DLMF 13.2.7: U(-m, b, 0) = (-1)^m (b)_m
                // (b)_m > m! overflows for m >= 171
                double r = 1.0;
                for (double k = 0.0; k < -a && k < 171.0 && std::isfinite(r); ++k) {
                    r *= b + k;
                }
                return odd ? -r : r;
            }
            // DLMF 13.2.16-18
            set_error("hyperu", SF_ERROR_SINGULAR, NULL);
            return odd ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
        }
        // DLMF 13.2.14-15 and 13.2.19-21
        return cephes::poch(1.0 - b + a, -a);
    }

    if (b == 1.0 && x < 1.0 && -0.25 < a && a < 0.3) {
        // DLMF 13.3.7. Fixes gh-15650
        return (x + 1.0 + 2.0 * a) * hypu(a + 1.0, 1.0, x) - std::pow(a + 1.0, 2) * hypu(a + 2.0, 1.0, x);
    }

    return hypu(a, b, x);
}

inline float hyperu(float a, float b, float x) {
    return static_cast<float>(hyperu(static_cast<double>(a), static_cast<double>(b), static_cast<double>(x)));
}

} // namespace xsf
