/* Translated from Cython into C++ by SciPy developers in 2024.
 *
 * Original authors: Pauli Virtanen, Eric Moore
 */

// Binomial coefficient

#pragma once

#include "config.h"

#include "cephes/beta.h"
#include "cephes/gamma.h"

namespace xsf {

XSF_HOST_DEVICE inline double binom(double n, double k) {
    double kx, nx, num, den, dk, sgn;

    if (n < 0) {
        nx = cxx::floor(n);
        if (n == nx) {
            // Undefined
            return cxx::numeric_limits<double>::quiet_NaN();
        }
    }

    kx = cxx::floor(k);
    if (k == kx && (cxx::abs(n) > 1E-8 || n == 0)) {
        /* Integer case: use multiplication formula for less rounding
         * error for cases where the result is an integer.
         *
         * This cannot be used for small nonzero n due to loss of
         * precision. */
        nx = cxx::floor(n);
        if (nx == n && kx > nx / 2 && nx > 0) {
            // Reduce kx by symmetry
            kx = nx - kx;
        }

        if (kx >= 0 && kx < 20) {
            num = 1.0;
            den = 1.0;
            for (int i = 1; i < 1 + static_cast<int>(kx); i++) {
                num *= i + n - kx;
                den *= i;
                if (cxx::abs(num) > 1E50) {
                    num /= den;
                    den = 1.0;
                }
            }
            return num / den;
        }
    }

    // general case
    if (n >= 1E10 * k && k > 0) {
        // avoid under/overflows intermediate results
        return cxx::exp(-cephes::lbeta(1 + n - k, 1 + k) - cxx::log(n + 1));
    }
    if (k > 1E8 * cxx::abs(n)) {
        // avoid loss of precision
        num = cephes::Gamma(1 + n) / cxx::abs(k) + cephes::Gamma(1 + n) * n / (2 * k * k); // + ...
        num /= M_PI * cxx::pow(cxx::abs(k), n);
        if (k > 0) {
            kx = cxx::floor(k);
            if (static_cast<int>(kx) == kx) {
                dk = k - kx;
                sgn = (static_cast<int>(kx) % 2 == 0) ? 1 : -1;
            } else {
                dk = k;
                sgn = 1;
            }
            return num * cxx::sin((dk - n) * M_PI) * sgn;
        }
        kx = cxx::floor(k);
        if (static_cast<int>(kx) == kx) {
            return 0;
        }
        return num * cxx::sin(k * M_PI);
    }
    return 1 / (n + 1) / cephes::beta(1 + n - k, 1 + k);
}

XSF_HOST_DEVICE inline float binom(float n, float k) { return binom(static_cast<double>(n), static_cast<double>(k)); }

} // namespace xsf
