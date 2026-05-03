/* Translated from Cython into C++ by SciPy developers in 2023.
 * Original header with Copyright information appears below.
 */

/* Implementation of the Lambert W function [1]. Based on MPMath
 *  Implementation [2], and documentation [3].
 *
 * Copyright: Yosef Meller, 2009
 * Author email: mellerf@netvision.net.il
 *
 * Distributed under the same license as SciPy
 *
 *
 * References:
 * [1] On the Lambert W function, Adv. Comp. Math. 5 (1996) 329-359,
 *     available online: https://web.archive.org/web/20230123211413/https://cs.uwaterloo.ca/research/tr/1993/03/W.pdf
 * [2] mpmath source code,
 https://github.com/mpmath/mpmath/blob/c5939823669e1bcce151d89261b802fe0d8978b4/mpmath/functions/functions.py#L435-L461
 * [3]
 https://web.archive.org/web/20230504171447/https://mpmath.org/doc/current/functions/powers.html#lambert-w-function
 *

 * TODO: use a series expansion when extremely close to the branch point
 * at `-1/e` and make sure that the proper branch is chosen there.
 */

#pragma once

#include "config.h"
#include "error.h"
#include "evalpoly.h"

namespace xsf {
constexpr double EXPN1 = 0.36787944117144232159553; // exp(-1)
constexpr double OMEGA = 0.56714329040978387299997; // W(1, 0)

namespace detail {
    template <typename T>
    XSF_HOST_DEVICE inline T lambertw_branchpt(T z) {
        // Series for W(z, 0) around the branch point; see 4.22 in [1].
        double coeffs[] = {-1.0 / 3.0, 1.0, -1.0};
        T p = std::sqrt(2.0 * (M_E * z + 1.0));

        return evalpoly(coeffs, 2, p);
    }

    template <typename T>
    XSF_HOST_DEVICE inline T lambertw_pade0(T z) {
        // (3, 2) Pade approximation for W(z, 0) around 0.
        double num[] = {12.85106382978723404255, 12.34042553191489361902, 1.0};
        double denom[] = {32.53191489361702127660, 14.34042553191489361702, 1.0};

        /* This only gets evaluated close to 0, so we don't need a more
         * careful algorithm that avoids overflow in the numerator for
         * large z. */
        return z * evalpoly(num, 2, z) / evalpoly(denom, 2, z);
    }

    template <typename T>
    XSF_HOST_DEVICE inline T lambertw_asy(T z, long k) {
        /* Compute the W function using the first two terms of the
         * asymptotic series. See 4.20 in [1].
         */
        T w;
        if constexpr(std::is_floating_point<T>::value) {
            w = std::log(z);
        } else {
            w = std::log(z) + 2.0 * M_PI * k * T(0, 1);
        }
        return w - std::log(w);
    }

    template <typename T>
    XSF_HOST_DEVICE inline T lambertw_impl(T z, long k, double tol) {
        double absz;
        double zreal, zimag, wreal;
        T w;
        T ew, wew, wewz, wn;

        if constexpr(std::is_floating_point<T>::value) {
            // real input
            zreal = z;
            zimag = 0;
        } else {
            // complex input
            zreal = z.real();
            zimag = z.imag();
        }

        // input is complex
        if (std::isnan(zreal) || std::isnan(zimag)) {
            return z;
        }
        if (zreal == std::numeric_limits<double>::infinity()) {
            if constexpr(std::is_floating_point<T>::value){
                return z;
            } else {
                return z + 2.0 * M_PI * k * std::complex<double>(0, 1);
            }
        }
        if (zreal == -std::numeric_limits<double>::infinity()) {
            if constexpr(std::is_floating_point<T>::value){
                return -z;
            } else {
                return -z + (2.0 * M_PI * k + M_PI) * std::complex<double>(0, 1);
            }
        }
        if (z == 0.0) {
            if (k == 0) {
                return z;
            }
            set_error("lambertw", SF_ERROR_SINGULAR, NULL);
            return -std::numeric_limits<double>::infinity();
        }
        if (z == 1.0 && k == 0) {
            // Split out this case because the asymptotic series blows up
            return OMEGA;
        }

        absz = std::abs(z);
        // Get an initial guess for Halley's method
        if (k == 0) {
            if (std::abs(z + EXPN1) < 0.3) {
                w = lambertw_branchpt(z);
            } else if (-1.0 < zreal && zreal < 1.5 && std::abs(zimag) < 1.0 &&
                    -2.5 * std::abs(zimag) - 0.2 < zreal) {
                /* Empirically determined decision boundary where the Pade
                * approximation is more accurate. */
                w = lambertw_pade0(z);
            } else {
                w = lambertw_asy(z, k);
            }
        } else if (k == -1) {
            if (absz <= EXPN1 && zimag == 0.0 && zreal < 0.0) {
                w = std::log(-zreal);
            } else {
                w = lambertw_asy(z, k);
            }
        } else {
            w = lambertw_asy(z, k);
        }

        if constexpr(std::is_floating_point<T>::value) {
            wreal = w;
        } else {
            wreal = w.real();
        }

        // Halley's method; see 5.9 in [1]
        if (wreal >= 0) {
            // Rearrange the formula to avoid overflow in exp
            for (int i = 0; i < 100; i++) {
                ew = std::exp(-w);
                wewz = w - z * ew;
                wn = w - wewz / (w + 1.0 - (w + 2.0) * wewz / (2.0 * w + 2.0));
                if (std::abs(wn - w) <= tol * std::abs(wn)) {
                    return wn;
                }
                w = wn;
            }
        } else {
            for (int i = 0; i < 100; i++) {
                ew = std::exp(w);
                wew = w * ew;
                wewz = wew - z;
                wn = w - wewz / (wew + ew - (w + 2.0) * wewz / (2.0 * w + 2.0));
                if (std::abs(wn - w) <= tol * std::abs(wn)) {
                    return wn;
                }
                w = wn;
            }
        }

        set_error("lambertw", SF_ERROR_SLOW, "iteration failed to converge: %g + %gj", zreal, zimag);
        if constexpr(std::is_floating_point<T>::value) {
            return std::numeric_limits<double>::quiet_NaN();
        } else {
            return std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }
    }

} // namespace detail

XSF_HOST_DEVICE inline double lambertw(double x, long k, double tol) {
    return detail::lambertw_impl(x, k, tol);
}

XSF_HOST_DEVICE inline float lambertw(float x, long k, float tol) {
    return static_cast<float>(detail::lambertw_impl(static_cast<double>(x), k, tol)
    );
}

XSF_HOST_DEVICE inline std::complex<double> lambertw(std::complex<double> z, long k, double tol) {
    return detail::lambertw_impl(z, k, tol);
}

XSF_HOST_DEVICE inline std::complex<float> lambertw(std::complex<float> z, long k, float tol) {
    return static_cast<std::complex<float>>(detail::lambertw_impl(static_cast<std::complex<double>>(z), k, static_cast<double>(tol))
    );
}

} // namespace xsf
