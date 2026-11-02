/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     ellik.c
 *
 *     Incomplete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double phi, m, y, ellik();
 *
 * y = ellik( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi | m) =   |    ------------------
 *                |                   2
 *              | |    sqrt( 1 - m sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points with m in [0, 1] and phi as indicated.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -10,10       200000      7.4e-16     1.0e-16
 *
 *
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
/* Copyright 2014, Eric W. Moore */

/*     Incomplete elliptic integral of first kind      */
#pragma once

#include "../config.h"
#include "../error.h"
#include "const.h"
#include "ellpk.h"

namespace xsf {
namespace cephes {

    namespace detail {

        /* To calculate legendre's incomplete elliptical integral of the first kind for
         * negative m, we use a power series in phi for small m*phi*phi, an asymptotic
         * series in m for large m*phi*phi* and the relation to Carlson's symmetric
         * integral of the first kind.
         *
         * F(phi, m) = sin(phi) * R_F(cos(phi)^2, 1 - m * sin(phi)^2, 1.0)
         *           = R_F(c-1, c-m, c)
         *
         * where c = csc(phi)^2. We use the second form of this for (approximately)
         * phi > 1/(sqrt(DBL_MAX) ~ 1e-154, where csc(phi)^2 overflows. Elsewhere we
         * use the first form, accounting for the smallness of phi.
         *
         * The algorithm used is described in Carlson, B. C. Numerical computation of
         * real or complex elliptic integrals. (1994) https://arxiv.org/abs/math/9409227
         * Most variable names reflect Carlson's usage.
         *
         * In this routine, we assume m < 0 and  0 > phi > pi/2.
         */
        XSF_HOST_DEVICE inline double ellik_neg_m(double phi, double m) {
            double x, y, z, x1, y1, z1, A0, A, Q, X, Y, Z, E2, E3, scale;
            int n = 0;
            double mpp = (m * phi) * phi;

            if (-mpp < 1e-6 && phi < -m) {
                return phi + (-mpp * phi * phi / 30.0 + 3.0 * mpp * mpp / 40.0 + mpp / 6.0) * phi;
            }

            if (-mpp > 4e7) {
                double sm = cxx::sqrt(-m);
                double sp = cxx::sin(phi);
                double cp = cxx::cos(phi);

                double a = cxx::log(4 * sp * sm / (1 + cp));
                double b = -(1 + cp / sp / sp - a) / 4 / m;
                return (a + b) / sm;
            }

            if (phi > 1e-153 && m > -1e305) {
                double s = cxx::sin(phi);
                double csc2 = 1.0 / (s * s);
                scale = 1.0;
                x = 1.0 / (cxx::tan(phi) * cxx::tan(phi));
                y = csc2 - m;
                z = csc2;
            } else {
                scale = phi;
                x = 1.0;
                y = 1 - m * scale * scale;
                z = 1.0;
            }

            if (x == y && x == z) {
                return scale / cxx::sqrt(x);
            }

            A0 = (x + y + z) / 3.0;
            A = A0;
            x1 = x;
            y1 = y;
            z1 = z;
            /* Carlson gives 1/pow(3*r, 1.0/6.0) for this constant. if r == eps,
             * it is ~338.38. */
            Q = 400.0 * cxx::fmax(cxx::abs(A0 - x), cxx::fmax(cxx::abs(A0 - y), cxx::abs(A0 - z)));

            while (Q > cxx::abs(A) && n <= 100) {
                double sx = cxx::sqrt(x1);
                double sy = cxx::sqrt(y1);
                double sz = cxx::sqrt(z1);
                double lam = sx * sy + sx * sz + sy * sz;
                x1 = (x1 + lam) / 4.0;
                y1 = (y1 + lam) / 4.0;
                z1 = (z1 + lam) / 4.0;
                A = (x1 + y1 + z1) / 3.0;
                n += 1;
                Q /= 4;
            }
            X = (A0 - x) / A / (1 << 2 * n);
            Y = (A0 - y) / A / (1 << 2 * n);
            Z = -(X + Y);

            E2 = X * Y - Z * Z;
            E3 = X * Y * Z;

            return scale * (1.0 - E2 / 10.0 + E3 / 14.0 + E2 * E2 / 24.0 - 3.0 * E2 * E3 / 44.0) / sqrt(A);
        }

    } // namespace detail

    XSF_HOST_DEVICE inline double ellik(double phi, double m) {
        double a, b, c, e, temp, t, K, denom, npio2;
        int d, mod, sign;

        if (cxx::isnan(phi) || cxx::isnan(m))
            return cxx::numeric_limits<double>::quiet_NaN();
        if (m > 1.0)
            return cxx::numeric_limits<double>::quiet_NaN();
        if (cxx::isinf(phi) || cxx::isinf(m)) {
            if (cxx::isinf(m) && cxx::isfinite(phi))
                return 0.0;
            else if (cxx::isinf(phi) && cxx::isfinite(m))
                return phi;
            else
                return cxx::numeric_limits<double>::quiet_NaN();
        }
        if (m == 0.0)
            return (phi);
        a = 1.0 - m;
        if (a == 0.0) {
            if (cxx::abs(phi) >= (double)M_PI_2) {
                set_error("ellik", SF_ERROR_SINGULAR, NULL);
                return (cxx::numeric_limits<double>::infinity());
            }
            /* DLMF 19.6.8, and 4.23.42 */
            return cxx::asinh(cxx::tan(phi));
        }
        npio2 = floor(phi / M_PI_2);
        if (cxx::fmod(cxx::abs(npio2), 2.0) == 1.0)
            npio2 += 1;
        if (npio2 != 0.0) {
            K = ellpk(a);
            phi = phi - npio2 * M_PI_2;
        } else
            K = 0.0;
        if (phi < 0.0) {
            phi = -phi;
            sign = -1;
        } else
            sign = 0;
        if (a > 1.0) {
            temp = detail::ellik_neg_m(phi, m);
            goto done;
        }
        b = cxx::sqrt(a);
        t = cxx::tan(phi);
        if (cxx::abs(t) > 10.0) {
            /* Transform the amplitude */
            e = 1.0 / (b * t);
            /* ... but avoid multiple recursions.  */
            if (cxx::abs(e) < 10.0) {
                e = cxx::atan(e);
                if (npio2 == 0)
                    K = ellpk(a);
                temp = K - ellik(e, m);
                goto done;
            }
        }
        a = 1.0;
        c = cxx::sqrt(m);
        d = 1;
        mod = 0;

        while (cxx::abs(c / a) > detail::MACHEP) {
            temp = b / a;
            phi = phi + atan(t * temp) + mod * M_PI;
            denom = 1.0 - temp * t * t;
            if (cxx::abs(denom) > 10 * detail::MACHEP) {
                t = t * (1.0 + temp) / denom;
                mod = (phi + M_PI_2) / M_PI;
            } else {
                t = cxx::tan(phi);
                mod = static_cast<int>(cxx::floor((phi - cxx::atan(t)) / M_PI));
            }
            c = (a - b) / 2.0;
            temp = cxx::sqrt(a * b);
            a = (a + b) / 2.0;
            b = temp;
            d += d;
        }

        temp = (cxx::atan(t) + mod * M_PI) / (d * a);

    done:
        if (sign < 0)
            temp = -temp;
        temp += npio2 * K;
        return (temp);
    }

} // namespace cephes
} // namespace xsf
