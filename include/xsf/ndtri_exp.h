// Translated from Cython to C++ by the xsf developers in 2026.

// Implementation of the inverse of the logarithm of the CDF of the standard
// normal distribution.

// Copyright: Albert Steppi

// Distributed under the same license as SciPy

// Implementation Overview

// The inverse of the CDF of the standard normal distribution is available
// in scipy through the Cephes Math Library where it is called ndtri.
// We call our implementation of the inverse of the log CDF ndtri_exp.
// For -2 <= y <= log(1 - exp(-2)),  ndtri_exp is computed as ndtri(exp(y)).

// For 0 < p < exp(-2), the cephes implementation of ndtri uses an approximation
// for ndtri(p) which is a function of z = sqrt(-2.0 * log(p)). Letting
// y = log(p), for y < -2, ndtri_exp uses this approximation in log(p) directly.
// This allows the implementation to achieve high precision for very small y,
// whereas ndtri(exp(y)) evaluates to infinity. This is because  exp(y) underflows
// for y < ~ -745.1.

// When p > 1 - exp(-2), the Cephes implementation of ndtri uses the symmetry
// of the normal distribution and calculates ndtri(p) as -ndtri(1 - p) allowing
// for the use of the same approximation. When y > log(1 - exp(-2)) this
// implementation calculates ndtri_exp as -ndtri(-expm1(y)).

// Accuracy

// Cephes provides the following relative error estimates for ndtri
//                      Relative error:
// arithmetic   domain        # trials      peak         rms
//    IEEE     0.125, 1        20000       7.2e-16     1.3e-16
//    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17

// When y < -2, ndtri_exp must have relative error at least as small as the
// Cephes implementation of ndtri for p < exp(-2). It relies on the same
// approximation but does not have to lose precision by passing from p to log(p)
// before applying the approximation.

// Relative error of ndtri for values of the argument p near 1 can be much higher
// than claimed by the above chart. For p near 1, symmetry is exploited to
// replace the calculation of ndtri(p) with -ndtri(1 - p). The inverse of the
// normal CDF increases so rapidly near the endpoints of [0, 1] that the loss
// of precision incurred by the subtraction 1 - p due to limitations in binary
// approximation can make a significant difference in the results. Using
// version 9.3.0 targeting x86_64-linux-gnu we've observed the following

//                                               Estimated Relative Error
//  ndtri(1e-8)      = -5.612001244174789        ''
// -ndtri(1 - 1e-8)  = -5.612001243305505        1.55e-10
//  ndtri(1e-16)     = -8.222082216130435        ''
// -ndtri(1 - 1e-16) = -8.209536151601387        0.0015

// If expm1 is correctly rounded for y in [log(1 - exp(-2), 0), then ndtri_exp(y)
// should have the same relative error as ndtri(p) for p > 1 - exp(-2). As seen
// above, this error may be higher than desired. IEEE-754 provides no guarantee on
// the accuracy of expm1 however, therefore accuracy of ndtri_exp in this range
// is platform dependent.

// The case

//     -2 <= y <= log(1 - exp(-2)) ~ -0.1454

// corresponds to

//      ~ 0.135 <= p <= ~ 0.865

// The derivative of ndtri is sqrt(2 * pi) * exp(ndtri(x)**2 / 2).
// It is ~4.597 at x ~ 0.135, decreases monotonically to sqrt(2 * pi) ~ 2.507
// at x = 0 and increases monotonically again to ~4.597 at x ~ 0.865.

// It can be checked that all higher derivatives follow a similar pattern.
// Their absolute value takes on a maximum (for this interval) at x ~ 0.135,
// decrease to a minimum at x = 0 and increases to the same maximum at x ~ 0.865.
// Derivatives of all orders are positive at x=log(1 - exp(-2)). Thus the worst
// possible loss of precision of ndtri(exp(x)) in the interval
// [0, log(1 - exp(-2))] due to error in calculating exp(x) must occur near
// x=log(1 - exp(-2)). By symmetry, the worst possible loss of precision in
// [-2, log(1 - exp(-2)] must occur near the endpoints. We may observe
// empirically that error at the endpoints due to exp is not substantial.
// Assuming that exp is accurate within +-ULP (unit of least precision),
// we observed a value of at most ~6.0474e-16 for

//     abs(ndtri(x + epsilon) - ndtri(x))

// if x is near exp(-2) or 1 - exp(-2) and epsilon is equal to the unit of least
// precision of x.

// (IEEE-754 provides no guarantee on the accuracy of exp, but for most
// compilers on most architectures an assumption of +-ULP should be
// reasonable.)

// The error here is on the order of the error in the Cephes implementation of
// ndtri itself, leading to an error profile that is still favorable.

#pragma once

#include "cephes/ndtri.h"
#include "config.h"

namespace xsf {

namespace detail {

    XSF_HOST_DEVICE inline double ndtri_exp_small_y(double y) {
        // Return the inverse of the logarithm of the normal CDF for very small y.

        // For p sufficiently small, the inverse of the CDF of the normal
        // distribution can be approximated to high precision as a rational function
        // in sqrt(-2.0 * log(p)).

        double x;
        // sqrt(-2 * y) is faster and has more precision but overflows when
        // y < -DBL_MAX * 0.5
        if (y >= -std::numeric_limits<double>::max() * 0.5) {
            x = std::sqrt(-2.0 * y);
        } else {
            x = M_SQRT2 * std::sqrt(-y);
        }

        const double x0 = x - std::log(x) / x;
        const double z = 1.0 / x;
        double x1;
        if (x < 8.0) {
            x1 = z * cephes::polevl(z, cephes::detail::ndtri_P1, 8) / cephes::p1evl(z, cephes::detail::ndtri_Q1, 8);
        } else {
            x1 = z * cephes::polevl(z, cephes::detail::ndtri_P2, 8) / cephes::p1evl(z, cephes::detail::ndtri_Q2, 8);
        }
        return x1 - x0;
    }

} // namespace detail

// Inverse of the logarithm of the normal CDF.
XSF_HOST_DEVICE inline double ndtri_exp(double y) {
    if (y < -std::numeric_limits<double>::max()) {
        return -std::numeric_limits<double>::infinity();
    }
    if (y < -2.0) {
        return detail::ndtri_exp_small_y(y);
    }
    if (y > -0.14541345786885906) { // log1p(-exp(-2))
        return -cephes::ndtri(-std::expm1(y));
    }
    return cephes::ndtri(std::exp(y));
}

XSF_HOST_DEVICE inline float ndtri_exp(float y) { return ndtri_exp(static_cast<double>(y)); }

} // namespace xsf
