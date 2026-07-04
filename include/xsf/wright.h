/**********************************************************************/
/* wrightomega is the simple routine for evaluating the wright omega  */
/* function.                                                          */
/*                                                                    */
/* Calling:                                                           */
/*    w = wrightomega(z)                                              */
/*                                                                    */
/* Input:                                                             */
/*   z  --  double complex                                            */
/*          Value to evaluate Wrightomega(z) at.                      */
/*                                                                    */
/* Output:                                                            */
/*   w  --  double complex                                            */
/*          Value of Wrightomega(z)                                   */
/*                                                                    */
/**********************************************************************/

/*
  Also published as ACM TOMS 917; relicensed as BSD by the author.

  Copyright (C) Piers Lawrence.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors
  may be used to endorse or promote products derived from this software without
  specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include "config.h"
#include "error.h"

namespace xsf {
namespace detail {

    constexpr double TWOITERTOL = std::numeric_limits<double>::epsilon();

    /* Computes fl(a+b) and err(a+b).  */
    XSF_HOST_DEVICE inline double two_sum(double a, double b, double *err) {
        volatile double s = a + b;
        volatile double c = s - a;
        volatile double d = b - c;
        volatile double e = s - c;
        *err = (a - e) + d;
        return s;
    }

    XSF_HOST_DEVICE inline double add_round_up(double a, double b) {
        double s, err;

        if (std::isnan(a) || std::isnan(b)) {
            return NAN;
        }

        s = two_sum(a, b, &err);
        if (err > 0) {
            /* fl(a + b) rounded down */
            return std::nextafter(s, INFINITY);
        } else {
            /* fl(a + b) rounded up or didn't round */
            return s;
        }
    }

    XSF_HOST_DEVICE inline double add_round_down(double a, double b) {
        double s, err;

        if (std::isnan(a) || std::isnan(b)) {
            return NAN;
        }

        s = two_sum(a, b, &err);
        if (err < 0) {
            return std::nextafter(s, -INFINITY);
        } else {
            return s;
        }
    }

} // namespace detail

XSF_HOST_DEVICE inline std::complex<double> wrightomega(std::complex<double> z) {
    constexpr std::complex<double> I(0.0, 1.0);
    double pi = M_PI, s = 1.0;
    double x, y, ympi, yppi, near;
    std::complex<double> e, r, pz, wp1, t, fac, w;

    /* extract real and imaginary parts of z */
    x = std::real(z);
    y = std::imag(z);

    /* compute if we are near the branch cuts */
    ympi = y - pi;
    yppi = y + pi;
    near = 0.1e-1;

    /* Test for floating point exceptions */
    /*****************************/
    /* NaN output for NaN input  */
    /*****************************/
    if (std::isnan(x) || std::isnan(y)) {
        return std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }
    /*********************************/
    /* Signed zeros between branches */
    /*********************************/
    else if (std::isinf(x) && (x < 0.0) && (-pi < y) && (y <= pi)) {
        if (std::fabs(y) <= pi / 2.0) {
            if (y >= 0) {
                return std::complex<double>(0.0, 0.0);
            } else {
                return std::complex<double>(0.0, -0.0);
            }
        } else {
            if (y >= 0) {
                return std::complex<double>(-0.0, 0.0);
            } else {
                return std::complex<double>(-0.0, -0.0);
            }
        }
    }
    /**************************/
    /* Asymptotic for large z */
    /**************************/
    else if (std::isinf(x) || std::isinf(y)) {
        return std::complex<double>(x, y);
    }

    /******************************************/
    /* Test If exactly on the singular points */
    /******************************************/
    if ((x == -1.0) && (std::fabs(y) == pi)) {
        return std::complex<double>(-1.0, 0.0);
    }

    /* Choose approximation based on region */
    /**********************************/
    /* Region 1: upper branch point   */
    /* Series about z=-1+Pi*I         */
    /**********************************/
    if ((-2.0 < x && x <= 1.0 && 1.0 < y && y < 2.0 * pi)) {
        pz = std::conj(std::sqrt(std::conj(2.0 * (z + 1.0 - I * pi))));
        w = -1.0 + (I + (1.0 / 3.0 + (-1.0 / 36.0 * I + (1.0 / 270.0 + 1.0 / 4320.0 * I * pz) * pz) * pz) * pz) * pz;
    }
    /**********************************/
    /* Region 2: lower branch point   */
    /* Series about z=-1-Pi*I         */
    /**********************************/
    else if ((-2.0 < x && x <= 1.0 && -2.0 * pi < y && y < -1.0)) {
        pz = std::conj(std::sqrt(std::conj(2.0 * (z + 1.0 + I * pi))));
        w = -1.0 + (-I + (1.0 / 3.0 + (1.0 / 36.0 * I + (1.0 / 270.0 - 1.0 / 4320.0 * I * pz) * pz) * pz) * pz) * pz;
    }
    /*********************************/
    /* Region 3: between branch cuts */
    /* Series: About -infinity       */
    /*********************************/
    else if (x <= -2.0 && -pi < y && y <= pi) {
        pz = std::exp(z);
        w = (1.0 + (-1.0 + (3.0 / 2.0 + (-8.0 / 3.0 + 125.0 / 24.0 * pz) * pz) * pz) * pz) * pz;
        if (w == 0.0) {
            set_error("wrightomega", SF_ERROR_UNDERFLOW, "underflow in exponential series");
            /* Skip the iterative scheme because it computes log(*w) */
            return w;
        }
    }
    /**********************/
    /* Region 4: Mushroom */
    /* Series about z=1   */
    /**********************/
    else if (
        ((-2.0 < x) && (x <= 1.0) && (-1.0 <= y) && (y <= 1.0)) ||
        ((-2.0 < x) && (x - 0.1e1) * (x - 0.1e1) + y * y <= pi * pi)
    ) {
        pz = z - 1.0;
        w = 1.0 / 2.0 + 1.0 / 2.0 * z +
            (1.0 / 16.0 + (-1.0 / 192.0 + (-1.0 / 3072.0 + 13.0 / 61440.0 * pz) * pz) * pz) * pz * pz;
    }
    /*************************/
    /* Region 5: Top wing    */
    /* Negative log series   */
    /*************************/
    else if (x <= -0.105e1 && pi < y && y - pi <= -0.75e0 * (x + 0.1e1)) {
        t = z - I * pi;
        pz = std::log(-t);
        w = t - pz;
        fac = pz / t;
        w += fac;
        fac /= t;
        w += fac * (0.5 * pz - 1.0);
        fac /= t;
        w += fac * (pz * pz / 3.0 - 3.0 * pz / 2.0 + 1.0);
        if (std::abs(z) > 1e50)
        /* Series is accurate and the iterative scheme could overflow */
        {
            return w;
        }
    }
    /***************************/
    /* Region 6: Bottom wing   */
    /* Negative log series     */
    /***************************/
    else if (x <= -0.105e1 && 0.75e0 * (x + 0.1e1) < y + pi && y + pi <= 0.0e0) {
        t = z + I * pi;
        pz = std::log(-t);
        w = t - pz;
        fac = pz / t;
        w += fac;
        fac /= t;
        w += fac * (0.5 * pz - 1.0);
        fac /= t;
        w += fac * (pz * pz / 3.0 - 3.0 * pz / 2.0 + 1.0);
        if (std::abs(z) > 1e50)
        /* Series is accurate and the iterative scheme could overflow */
        {
            return w;
        }
    }
    /************************************/
    /* Region 7: Everywhere else        */
    /* Series solution about infinity   */
    /************************************/
    else {
        pz = std::log(z);
        w = z - pz;
        fac = pz / z;
        w += fac;
        fac /= z;
        w += fac * (0.5 * pz - 1.0);
        fac /= z;
        w += fac * (pz * pz / 3.0 - 3.0 * pz / 2.0 + 1.0);
        if (std::abs(z) > 1e50)
        /* Series is accurate and the iterative scheme could overflow */
        {
            return w;
        }
    }

    /**********************************/
    /* Regularize if near branch cuts */
    /**********************************/
    if (x <= -0.1e1 + near && (std::fabs(ympi) <= near || std::fabs(yppi) <= near)) {
        s = -1.0;
        if (std::fabs(ympi) <= near) {
            /* Recompute ympi with directed rounding */
            ympi = detail::add_round_up(y, -pi);

            if (ympi <= 0.0) {
                ympi = detail::add_round_down(y, -pi);
            }

            z = x + I * ympi;
        } else {
            /* Recompute yppi with directed rounding */
            yppi = detail::add_round_up(y, pi);

            if (yppi <= 0.0) {
                yppi = detail::add_round_down(y, pi);
            }

            z = x + I * yppi;
        }
    }

    /*****************/
    /* Iteration one */
    /*****************/
    w = s * w;
    r = z - s * w - std::log(w);
    wp1 = s * w + 1.0;
    e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r) / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
    w = w * (1.0 + e);

    /*****************/
    /* Iteration two */
    /*****************/
    if (std::abs((2.0 * w * w - 8.0 * w - 1.0) * std::pow(std::abs(r), 4.0)) >=
        detail::TWOITERTOL * 72.0 * std::pow(std::abs(wp1), 6.0)) {
        r = z - s * w - std::log(w);
        wp1 = s * w + 1.0;
        e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r) / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
        w = w * (1.0 + e);
    }

    /***********************/
    /* Undo regularization */
    /***********************/
    w = s * w;

    return w;
}

XSF_HOST_DEVICE inline std::complex<float> wrightomega(std::complex<float> z) {
    return static_cast<std::complex<float>>(wrightomega(static_cast<std::complex<double>>(z)));
}

XSF_HOST_DEVICE inline double wrightomega(double x) {
    double w, wp1, e, r;

    /* NaN output for NaN input  */
    if (std::isnan(x)) {
        return x;
    }

    /* Positive infinity is asymptotically x */
    /* Negative infinity is zero */
    if (std::isinf(x)) {
        if (x > 0.0) {
            return x;
        } else {
            return 0.0;
        }
    }

    if (x < -50.0) {
        /*
         * Skip the iterative scheme because  exp(x) is already
         * accurate to double precision.
         */
        w = std::exp(x);
        if (w == 0.0) {
            set_error("wrightomega", SF_ERROR_UNDERFLOW, "underflow in exponential series");
        }
        return w;
    }
    if (x > 1e20) {
        /*
         * Skip the iterative scheme because the result is just x to
         * double precision
         */
        return x;
    }

    /* Split into three distinct intervals (-inf,-2), [-2,1), [1,inf) */
    if (x < -2.0) {
        /* exponential is approx < 1.3e-1 accurate */
        w = std::exp(x);
    } else if (x < 1) {
        /* on [-2,1) approx < 1.5e-1 accurate */
        w = std::exp(2.0 * (x - 1.0) / 3.0);
    } else {
        /* infinite series with 2 terms approx <1.7e-1 accurate */
        w = std::log(x);
        w = x - w + w / x;
    }

    /* Iteration one of Fritsch, Shafer, and Crowley (FSC) iteration */
    r = x - w - std::log(w);
    wp1 = w + 1.0;
    e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r) / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
    w = w * (1.0 + e);

    /* Iteration two (if needed based on the condition number) */
    if (std::fabs((2.0 * w * w - 8.0 * w - 1.0) * std::pow(std::fabs(r), 4.0)) >=
        detail::TWOITERTOL * 72.0 * std::pow(std::fabs(wp1), 6.0)) {
        r = x - w - std::log(w);
        wp1 = w + 1.0;
        e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r) / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
        w = w * (1.0 + e);
    }

    return w;
}

XSF_HOST_DEVICE inline float wrightomega(float x) { return static_cast<float>(wrightomega(static_cast<double>(x))); }
} // namespace xsf
