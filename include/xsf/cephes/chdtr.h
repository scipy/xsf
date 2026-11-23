/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     chdtr.c
 *
 *     Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double df, x, y, chdtr();
 *
 * y = chdtr( df, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the left hand tail (from 0 to x)
 * of the Chi square probability density function with
 * v degrees of freedom.
 *
 *
 *                                  inf.
 *                                    -
 *                        1          | |  v/2-1  -t/2
 *  P( x | v )   =   -----------     |   t      e     dt
 *                    v/2  -       | |
 *                   2    | (v/2)   -
 *                                   x
 *
 * where x is the Chi-square variable.
 *
 * The incomplete Gamma integral is used, according to the
 * formula
 *
 *     y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).
 *
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 * See igam().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * chdtr domain   x < 0 or v < 1        0.0
 */
/*							chdtrc()
 *
 *	Complemented Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double v, x, y, chdtrc();
 *
 * y = chdtrc( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the right hand tail (from x to
 * infinity) of the Chi square probability density function
 * with v degrees of freedom:
 *
 *
 *                                  inf.
 *                                    -
 *                        1          | |  v/2-1  -t/2
 *  P( x | v )   =   -----------     |   t      e     dt
 *                    v/2  -       | |
 *                   2    | (v/2)   -
 *                                   x
 *
 * where x is the Chi-square variable.
 *
 * The incomplete Gamma integral is used, according to the
 * formula
 *
 *	y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).
 *
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 * See igamc().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * chdtrc domain  x < 0 or v < 1        0.0
 */
/*							chdtri()
 *
 *	Inverse of complemented Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double df, x, y, chdtri();
 *
 * x = chdtri( df, y );
 *
 *
 *
 *
 * DESCRIPTION:
 *
 * Finds the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given cumulative probability y.
 *
 * This is accomplished using the inverse Gamma integral
 * function and the relation
 *
 *    x/2 = igamci( df/2, y );
 *
 *
 *
 *
 * ACCURACY:
 *
 * See igami.c.
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * chdtri domain   y < 0 or y > 1        0.0
 *                     v < 1
 *
 */

/*                                                             chdtr() */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"
#include "../error.h"

#include "igam.h"
#include "igami.h"

namespace xsf {
namespace cephes {

    XSF_HOST_DEVICE inline double chdtrc(double df, double x) {
        if (x < 0.0) {
            set_error("chdtr", SF_ERROR_DOMAIN, NULL);
            return (cxx::numeric_limits<double>::quiet_NaN());
        }
        if (x == 0.0 && !cxx::isnan(df) && df > 0.0) {
            return 1.0; // also if df / 2 underflows, as igamc(0, 0) is NaN
        }
        return (igamc(df / 2.0, x / 2.0));
    }

    XSF_HOST_DEVICE inline double chdtr(double df, double x) {

        if ((x < 0.0)) { /* || (df < 1.0) ) */
            set_error("chdtr", SF_ERROR_DOMAIN, NULL);
            return (cxx::numeric_limits<double>::quiet_NaN());
        }
        if (x < 2 * cxx::numeric_limits<double>::min() && !cxx::isnan(df) && df > 0.0) {
            if (x == 0.0) {
                return 0.0; // also if df / 2 underflows, as igam(0, 0) is NaN
            }
            // x / 2 would be subnormal; leading term of DLMF 8.7.1, (x / 2)^a / Gamma(a + 1)
            double a = df / 2.0;
            // split x^a so that the product rounds only once when it is subnormal
            double p = cxx::pow(x, a / 2.0);
            return p * (p * (cxx::pow(2.0, -a) / Gamma(a + 1.0)));
        }
        return (igam(df / 2.0, x / 2.0));
    }

    XSF_HOST_DEVICE inline double chdtri(double df, double y) {
        double x;

        if ((y < 0.0) || (y > 1.0)) { /* || (df < 1.0) ) */
            set_error("chdtri", SF_ERROR_DOMAIN, NULL);
            return (cxx::numeric_limits<double>::quiet_NaN());
        }

        x = igamci(0.5 * df, y);
        return (2.0 * x);
    }

} // namespace cephes
} // namespace xsf
