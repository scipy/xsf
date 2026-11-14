#pragma once

#include <cmath>

#include "cephes/dd_real.h"
#include "config.h"
#include "error.h"

namespace xsf {

template <typename T>
T expit(T x) {
    return 1 / (1 + std::exp(-x));
};

inline double exprel(double x) {
    if (std::abs(x) < std::numeric_limits<double>::epsilon()) {
        return 1;
    }

    if (x > 717) { // near log(DBL_MAX)
        return std::numeric_limits<double>::infinity();
    }

    return std::expm1(x) / x;
}

inline float exprel(float x) { return exprel(static_cast<double>(x)); }

template <typename T>
T logit(T x) {
    // the standard formula is log(x/(1 - x)), but this expression
    // loses precision near x=0.5, as does log(x) - log1p(-x)
    // we use the standard formula away from p=0.5, and use
    // log1p(2*(x - 0.5)) - log1p(-2*(x - 0.5)) around p=0.5, which
    // provides very good precision in this interval
    if (x < 0.3 || x > 0.65) {
        return std::log(x / (1 - x));
    } else {
        T s = 2 * (x - 0.5);
        return std::log1p(s) - std::log1p(-s);
    }
};

//
// the logistic sigmoid function 'expit' is
//
//     S(x) = 1/(1 + exp(-x))     = exp(x)/(exp(x) + 1)
//
// so
//
// log S(x) = -log(1 + exp(-x))   = x - log(exp(x) + 1)
//          = -log1p(exp(-x))     = x - log1p(exp(x))
//
// by using -log1p(exp(-x)) for x >= 0 and x - log1p(exp(x))
// for x < 0, we extend the range of x values for which we
// obtain accurate results (compared to the naive implementation
// log(expit(x)))
//
template <typename T>
T log_expit(T x) {
    if (x < 0) {
        return x - std::log1p(std::exp(x));
    }

    return -std::log1p(std::exp(-x));
};

/* compute log(1 - exp(x)) */
template <typename T>
T log1mexp(T x) {
    if (x > 0) {
        set_error("_log1mexp", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<T>::quiet_NaN();
    }
    if (x == 0) {
        set_error("_log1mexp", SF_ERROR_SINGULAR, NULL);
        return -std::numeric_limits<T>::infinity();
    }
    if (x < -1) {
        return std::log1p(-std::exp(x));
    }
    return std::log(-std::expm1(x));
}

/* compute log(exp(a) - exp(b)) for a >= b */
template <typename T>
T logdiffexp(T a, T b) {
    const T inf = std::numeric_limits<T>::infinity();

    if (a < b || (a == inf && b == inf)) {
        set_error("logdiffexp", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<T>::quiet_NaN();
    }
    if (a == b) {
        set_error("logdiffexp", SF_ERROR_SINGULAR, NULL);
        return -inf;
    }

    const T x = b - a;
    const T result = a + log1mexp(x);

    // limit the fallback to cancellation near a zero result
    if (!(a > 0 && std::abs(result) < T(0.5) * std::min(a, T(1)))) {
        return result;
    }

    using cephes::detail::double_double;

    const double aa = static_cast<double>(a);
    const double bb = static_cast<double>(b);

    if (x <= -700) {
        // scale the tail subtraction before rounding to subnormal values
        // the omitted quadratic term is below the double subnormal range
        const double_double d = double_double(bb) - double_double(aa);
        const double_double scaled_a = cephes::detail::ldexp(double_double(aa), 512);
        const double_double scaled_exp = exp(d + cephes::detail::LOG2 * 512.0);
        const double_double residual = scaled_a - scaled_exp;

        return static_cast<T>(std::ldexp(static_cast<double>(residual), -512));
    }

    double_double expm1_a;

    if (aa < 0.5) {
        // evaluate expm1(a) without subtracting nearly equal values
        double_double term(aa);
        expm1_a = term;

        for (int k = 2; k <= 40; ++k) {
            term = term * double_double(aa) / double_double(k);
            expm1_a = expm1_a + term;

            if (std::abs(static_cast<double>(term)) <= std::abs(static_cast<double>(expm1_a)) * 1e-32) {
                break;
            }
        }
    } else {
        expm1_a = exp(double_double(aa)) - double_double(1.0);
    }

    // compute exp(a) - exp(b) - 1 with extra precision
    // log1p then preserves the small residual near a zero result
    const double_double residual = expm1_a - exp(double_double(bb));

    return static_cast<T>(std::log1p(static_cast<double>(residual)));
}

} // namespace xsf
