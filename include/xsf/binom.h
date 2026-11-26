/* Translated from Cython into C++ by SciPy developers in 2024.
 *
 * Original authors: Pauli Virtanen, Eric Moore
 */

// Binomial coefficient

#pragma once

#include "config.h"

#include "cephes/beta.h"
#include "cephes/gamma.h"
#include "cephes/poch.h"
#include "cephes/trig.h"

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
        // avoid loss of precision (DLMF 5.11.13)
        num = cephes::Gamma(1 + n) / cxx::abs(k) * (1 + n * (n + 1) / (2 * k)); // + ...
        // split k^n, which can overflow when the result doesn't
        den = cxx::pow(cxx::abs(k), n / 2);
        num = num / (M_PI * den) / den;
        if (k > 0) {
            kx = cxx::floor(k);
            dk = k - kx;
            sgn = (cxx::fmod(kx, 2.0) == 0) ? 1 : -1;
            return num * cephes::sinpi(dk - n) * sgn;
        }
        kx = cxx::floor(k);
        if (static_cast<int>(kx) == kx) {
            return 0;
        }
        return num * cxx::sin(k * M_PI);
    }
    if (cxx::abs(k) > 2E4 && cxx::abs(n) < 170 && cxx::isfinite(k)) {
        // beta loses ~|k| eps here, so use the reflection formula
        // binom(n, k) = Gamma(1 + n) sin(pi (n - k + 1)) / (pi poch(a, n + 1)) with a = k - n,
        // or with a = -k and sin(-pi k) for k < 0 (DLMF 5.5.3)
        double a, s;
        if (k > 0) {
            a = k - n;
            // sin(pi (n - k + 1)) = (-1)^kx sin(pi ((k - kx) - n)), where k - kx is exact
            kx = cxx::round(a);
            sgn = (cxx::fmod(kx, 2.0) == 0) ? 1 : -1;
            s = sgn * cephes::sinpi((k - kx) - n);
        } else {
            a = -k;
            s = cephes::sinpi(-k);
        }
        if (s == 0) {
            return 0;
        }
        // poch(a, n + 1) = poch(a, n) (a + n), avoiding the rounding of n + 1
        num = cephes::poch(a, n);
        den = num * (a + n);
        if (cxx::abs(num) >= cxx::numeric_limits<double>::min() && cxx::isfinite(den)) {
            return cephes::Gamma(1 + n) / den * (s / M_PI);
        }
        // split poch(a, n) if it over- or underflows
        num = cephes::poch(a, n / 2);
        den = cephes::poch(a + n / 2, n / 2) * (a + n);
        return cephes::Gamma(1 + n) / num * (s / M_PI) / den;
    }
    return 1 / (n + 1) / cephes::beta(1 + n - k, 1 + k);
}

XSF_HOST_DEVICE inline float binom(float n, float k) { return binom(static_cast<double>(n), static_cast<double>(k)); }

} // namespace xsf
