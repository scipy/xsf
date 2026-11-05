/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     beta.c
 *
 *     Beta function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, y, beta();
 *
 * y = beta( a, b );
 *
 *
 *
 * DESCRIPTION:
 *
 *                   -     -
 *                  | (a) | (b)
 * beta( a, b )  =  -----------.
 *                     -
 *                    | (a+b)
 *
 * For large arguments the logarithm of the function is
 * evaluated using lgam(), then exponentiated.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,30       30000       8.1e-14     1.1e-14
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * beta overflow    log(beta) > MAXLOG       0.0
 *                  a or b <0 integer        0.0
 *
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"
#include "const.h"
#include "gamma.h"
#include "rgamma.h"

namespace xsf {
namespace cephes {

    XSF_HOST_DEVICE double beta(double, double);
    XSF_HOST_DEVICE double lbeta(double, double);

    namespace detail {
        /*
         * Stable Stirling difference for ln Gamma(a) - ln Gamma(a + b)
         * when a >= b > 0 and a is large (a > MAXGAM).
         *
         * Evaluating ln Gamma(a) - ln Gamma(a + b) directly cancels significant
         * digits when a is large, yielding relative error ~1e-9.
         * By expanding the difference analytically using Stirling's series:
         *   ln Gamma(a) - ln Gamma(a + b) = -b * ln(a) + Delta(a, b)
         * where:
         *   Delta(a, b) = b - (a + b - 1/2) * ln(1 + b/a) + sum_k [B_{2k}/(2k(2k-1))] * [1/a^{2k-1} - 1/(a+b)^{2k-1}]
         * For small x = b/a, b - a * ln(1 + x) = b * (x/2 - x^2/3 + x^3/4 - ...)
         * which avoids any cancellation.
         */
        XSF_HOST_DEVICE inline double lgam_diff_stirling(double a, double b) {
            double x = b / a;
            double term_log = -b * cxx::log(a);
            double cancellation;

            if (x < 1e-4) {
                cancellation = b * (x / 2.0 - (x * x) / 3.0 + (x * x * x) / 4.0 -
                                   (x * x * x * x) / 5.0 + (x * x * x * x * x) / 6.0);
            } else {
                cancellation = b - a * cxx::log1p(x);
            }

            double term_mid = cancellation - (b - 0.5) * cxx::log1p(x);

            /* First Stirling correction: 1/(12*a) - 1/(12*(a+b)) = (b/a) / (12 * (a + b)) */
            double corr1 = (a < 1e150) ? (x / (12.0 * (a + b))) : 0.0;

            /* Second Stirling correction: -1/(360*a^3) + 1/(360*(a+b)^3) */
            double corr2 = 0.0;
            if (a < 1e4) {
                double u = 1.0 + x;
                corr2 = -((u * u * u) - 1.0) / (360.0 * (a * a * a) * (u * u * u));
            }

            return term_log + term_mid + corr1 + corr2;
        }

        /*
         * Asymptotic expansion for ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
         */
        XSF_HOST_DEVICE inline double lbeta_asymp(double a, double b, int *sgn) {
            double r = lgam_sgn(b, sgn);
            r -= b * cxx::log(a);

            /*
             * Evaluate the correction terms through the small ratios b/a and (1-b)/a
             * instead of the products b*(1-b), b*b*(1-b)*(1-b) and a*a*a, which overflow
             * for |b| > ~1e77 or a > ~5.6e102 and turned the result into NaN (inf/inf)
             * for e.g. lbeta(1e100, 1e107) although log|B| is an ordinary finite number.
             */
            double t = b / a;
            double u = (1 - b) / a;

            r += t * (1 - b) / 2;
            r += t * u * (1 - 2 * b) / 12;
            r += -t * t * u * (1 - b) / 12;
            return r;
        }

        /*
         * Special case for a negative integer argument
         */

        XSF_HOST_DEVICE inline double beta_negint(int a, double b) {
            int sgn;
            if (b == static_cast<int>(b) && 1 - a - b > 0) {
                sgn = (static_cast<int>(b) % 2 == 0) ? 1 : -1;
                return sgn * xsf::cephes::beta(1 - a - b, b);
            } else {
                set_error("lbeta", SF_ERROR_OVERFLOW, NULL);
                return cxx::numeric_limits<double>::infinity();
            }
        }

        XSF_HOST_DEVICE inline double lbeta_negint(int a, double b) {
            double r;
            if (b == static_cast<int>(b) && 1 - a - b > 0) {
                r = xsf::cephes::lbeta(1 - a - b, b);
                return r;
            } else {
                set_error("lbeta", SF_ERROR_OVERFLOW, NULL);
                return cxx::numeric_limits<double>::infinity();
            }
        }
    } // namespace detail

    XSF_HOST_DEVICE inline double beta(double a, double b) {
        double y;
        int sign = 1;

        if (a <= 0.0) {
            if (a == cxx::floor(a)) {
                if (a == static_cast<int>(a)) {
                    return detail::beta_negint(static_cast<int>(a), b);
                } else {
                    goto overflow;
                }
            }
        }

        if (b <= 0.0) {
            if (b == cxx::floor(b)) {
                if (b == static_cast<int>(b)) {
                    return detail::beta_negint(static_cast<int>(b), a);
                } else {
                    goto overflow;
                }
            }
        }

        if (cxx::abs(a) < cxx::abs(b)) {
            y = a;
            a = b;
            b = y;
        }

        if (cxx::abs(a) > detail::beta_ASYMP_FACTOR * cxx::abs(b) && a > detail::beta_ASYMP_FACTOR) {
            /* Avoid loss of precision in lgam(a + b) - lgam(a) */
            y = detail::lbeta_asymp(a, b, &sign);
            return sign * cxx::exp(y);
        }

        y = a + b;
        if (cxx::abs(y) > detail::MAXGAM || cxx::abs(a) > detail::MAXGAM || cxx::abs(b) > detail::MAXGAM) {
            int sgngam;
            y = detail::lgam_sgn(y, &sgngam);
            sign *= sgngam; /* keep track of the sign */
            y = detail::lgam_sgn(b, &sgngam) - y;
            sign *= sgngam;
            y = detail::lgam_sgn(a, &sgngam) + y;
            sign *= sgngam;
            if (y > detail::MAXLOG) {
                goto overflow;
            }
            return (sign * cxx::exp(y));
        }

        y = rgamma(y);
        a = Gamma(a);
        b = Gamma(b);
        if (cxx::isinf(y)) {
            goto overflow;
        }

        if (cxx::abs(cxx::abs(a * y) - 1.0) > cxx::abs(cxx::abs(b * y) - 1.0)) {
            y = b * y;
            y *= a;
        } else {
            y = a * y;
            y *= b;
        }

        return (y);

    overflow:
        set_error("beta", SF_ERROR_OVERFLOW, NULL);
        return (sign * cxx::numeric_limits<double>::infinity());
    }

    /* Natural log of |beta|. */

    XSF_HOST_DEVICE inline double lbeta(double a, double b) {
        double y;
        int sign;

        sign = 1;

        if (a <= 0.0) {
            if (a == cxx::floor(a)) {
                if (a == static_cast<int>(a)) {
                    return detail::lbeta_negint(static_cast<int>(a), b);
                } else {
                    goto over;
                }
            }
        }

        if (b <= 0.0) {
            if (b == cxx::floor(b)) {
                if (b == static_cast<int>(b)) {
                    return detail::lbeta_negint(static_cast<int>(b), a);
                } else {
                    goto over;
                }
            }
        }

        if (cxx::abs(a) < cxx::abs(b)) {
            y = a;
            a = b;
            b = y;
        }

        if (a > detail::MAXGAM && b > 0.0) {
            /*
             * For large positive arguments, avoid catastrophic cancellation in
             * lgam(a) + lgam(b) - lgam(a + b) by using the stable Stirling difference.
             */
            int sgngam;
            double lgb = detail::lgam_sgn(b, &sgngam);
            return lgb + detail::lgam_diff_stirling(a, b);
        }

        y = a + b;
        if (cxx::abs(y) > detail::MAXGAM || cxx::abs(a) > detail::MAXGAM || cxx::abs(b) > detail::MAXGAM) {
            int sgngam;
            y = detail::lgam_sgn(y, &sgngam);
            sign *= sgngam; /* keep track of the sign */
            y = detail::lgam_sgn(b, &sgngam) - y;
            sign *= sgngam;
            y = detail::lgam_sgn(a, &sgngam) + y;
            sign *= sgngam;
            return (y);
        }

        y = rgamma(y);
        a = Gamma(a);
        b = Gamma(b);
        if (cxx::isinf(y)) {
        over:
            set_error("lbeta", SF_ERROR_OVERFLOW, NULL);
            return (sign * cxx::numeric_limits<double>::infinity());
        }

        if (cxx::abs(cxx::abs(a * y) - 1.0) > cxx::abs(cxx::abs(b * y) - 1.0)) {
            y = b * y;
            y *= a;
        } else {
            y = a * y;
            y *= b;
        }

        if (y < 0) {
            y = -y;
        }

        return (cxx::log(y));
    }
} // namespace cephes
} // namespace xsf
