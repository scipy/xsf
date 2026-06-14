/*
 * Implement Spence's function, a.k.a. the dilogarithm, for complex
 * arguments. Note that our definition differs from that in the sources
 * by the mapping z -> 1 - z.
 *
 * Sources
 * [1] Zagier, "The Dilogarithm Function"
 * [2] functions.wolfram.com
 * [3] Ginsberg, Zaborowski, "The Dilogarithm Function of a Real Argument"
 *
 * Author: Josh Wilson
 */

#pragma once

#include "config.h"
#include "zlog1.h"

namespace xsf {
namespace detail {
    // Relative tolerance for the series
    constexpr double TOL = 2.220446092504131e-16;
    constexpr double PISQ_6 = 1.6449340668482264365;

    XSF_HOST_DEVICE inline std::complex<double> cspence_series0(std::complex<double> z) {
        // A series centered at z = 0; see http://functions.wolfram.com/10.07.06.0005.02
        std::complex<double> zfac = 1;
        std::complex<double> sum1 = 0;
        std::complex<double> sum2 = 0;
        std::complex<double> term1, term2;

        if (z == 0.) {
            return PISQ_6;
        }

        for (int n = 1; n < 500; n++) {
            zfac *= z;
            term1 = zfac / static_cast<double>(n * n);
            sum1 += term1;
            term2 = zfac / static_cast<double>(n);
            sum2 += term2;
            if (std::abs(term1) <= TOL * std::abs(sum1) && std::abs(term2) <= TOL * std::abs(sum2)) {
                break;
            }
        }

        return (PISQ_6 - sum1 + detail::zlog1(z) * sum2);
    }

    XSF_HOST_DEVICE inline std::complex<double> cspence_series1(std::complex<double> z) {
        /*
         * A series centered at z = 1 which enjoys faster convergence than
         * the Taylor series. See [3]. The number of terms used comes from
         * bounding the absolute tolerance at the edge of the radius of
         * convergence where the sum is O(1).
         */
        std::complex<double> zfac = 1;
        std::complex<double> res = 0;
        std::complex<double> term, zz;

        if (z == 1.) {
            return 0;
        }

        z = 1. - z;
        zz = z * z;
        for (int n = 1; n < 500; n++) {
            zfac *= z;
            // Do the divisions one at a time to guard against overflow
            double dn = static_cast<double>(n);
            term = ((zfac / (dn * dn)) / ((dn + 1) * (dn + 1))) / ((dn + 2) * (dn + 2));
            res += term;
            if (std::abs(term) <= TOL * std::abs(res)) {
                break;
            }
        }
        res *= 4. * zz;
        res += 4. * z + 5.75 * zz + 3. * (1. - zz) * detail::zlog1(1. - z);
        res /= 1. + 4. * z + zz;
        return res;
    }

} // namespace detail

XSF_HOST_DEVICE inline std::complex<double> cspence(std::complex<double> z) {
    /*
     * Compute Spence's function for complex arguments. The strategy is:
     * - If z is close to 0, use a series centered at 0.
     * - If z is far away from 1, use the reflection formula
     *
     * spence(z) = -spence(z/(z - 1)) - pi**2/6 - ln(z - 1)**2/2
     *
     * to move close to 1. See [1].
     * - If z is close to 1, use a series centered at 1.
     */
    if (std::isnan(z.real()) || std::isnan(z.imag())) {
        return std::complex<double>{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    } else if (std::abs(z) < 0.5) {
        // This step isn't necessary, but this series converges faster.
        return detail::cspence_series0(z);
    } else if (std::abs(z - 1.) > 1.) {
        return -detail::cspence_series1(z / (z - 1.)) - detail::PISQ_6 -
               0.5 * detail::zlog1(z - 1.) * detail::zlog1(z - 1.);
    } else {
        return detail::cspence_series1(z);
    }
}

XSF_HOST_DEVICE inline std::complex<float> cspence(std::complex<float> z) {
    return static_cast<std::complex<float>>(cspence(static_cast<std::complex<double>>(z)));
}

} // namespace xsf
