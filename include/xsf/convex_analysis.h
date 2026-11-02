// Translated from Cython to C++ by the xsf developers in 2026.
// Original implementation by argriffing.
#pragma once

#include "config.h"

namespace xsf {

// Elementwise function for computing entropy.
XSF_HOST_DEVICE inline double entr(double x) {
    if (cxx::isnan(x)) {
        return x;
    } else if (x > 0) {
        return -x * cxx::log(x);
    } else if (x == 0) {
        return 0;
    } else {
        return -cxx::numeric_limits<double>::infinity();
    }
}

XSF_HOST_DEVICE inline float entr(float x) { return entr(static_cast<double>(x)); }

// Elementwise function for computing Kullback-Leibler divergence.
XSF_HOST_DEVICE inline double kl_div(double x, double y) {
    if (cxx::isnan(x) || cxx::isnan(y)) {
        return cxx::numeric_limits<double>::quiet_NaN();
    } else if (x > 0 && y > 0) {
        return x * cxx::log(x / y) - x + y;
    } else if (x == 0 && y >= 0) {
        return y;
    } else {
        return cxx::numeric_limits<double>::infinity();
    }
}

XSF_HOST_DEVICE inline float kl_div(float x, float y) { return kl_div(static_cast<double>(x), static_cast<double>(y)); }

// Elementwise function for computing relative entropy.
XSF_HOST_DEVICE inline double rel_entr(double x, double y) {
    if (cxx::isnan(x) || cxx::isnan(y)) {
        return cxx::numeric_limits<double>::quiet_NaN();
    }
    if (x <= 0 || y <= 0) {
        if (x == 0 && y >= 0) {
            return 0;
        }
        return cxx::numeric_limits<double>::infinity();
    }

    double ratio = x / y;
    if (0.5 < ratio && ratio < 2) {
        // When x and y are close, this is more accurate
        return x * cxx::log1p((x - y) / y);
    }
    if (cxx::numeric_limits<double>::min() < ratio && ratio < cxx::numeric_limits<double>::infinity()) {
        // There are no underflow/overflow issues
        return x * cxx::log(ratio);
    }
    // x and y are so far apart that taking x / y
    // results in either an underflow, overflow,
    // or subnormal number. Do the logarithm first
    return x * (cxx::log(x) - cxx::log(y));
}

XSF_HOST_DEVICE inline float rel_entr(float x, float y) {
    return rel_entr(static_cast<double>(x), static_cast<double>(y));
}

// Huber loss function.
XSF_HOST_DEVICE inline double huber(double delta, double r) {
    if (delta < 0) {
        return cxx::numeric_limits<double>::infinity();
    } else if (cxx::fabs(r) <= delta) {
        return 0.5 * r * r;
    } else {
        return delta * (cxx::fabs(r) - 0.5 * delta);
    }
}

XSF_HOST_DEVICE inline float huber(float delta, float r) {
    return huber(static_cast<double>(delta), static_cast<double>(r));
}

// Pseudo-Huber loss function.
XSF_HOST_DEVICE inline double pseudo_huber(double delta, double r) {
    if (delta < 0) {
        return cxx::numeric_limits<double>::infinity();
    } else if (delta == 0 || r == 0) {
        return 0;
    } else {
        double u = delta;
        double v = r / delta;
        // The formula is u*u*(sqrt(1 + v*v) - 1), but to maintain
        // precision with small v, we use
        //   sqrt(1 + v*v) - 1  =  exp(0.5*log(1 + v*v)) - 1
        //                      =  expm1(0.5*log1p(v*v))
        return u * u * cxx::expm1(0.5 * cxx::log1p(v * v));
    }
}

XSF_HOST_DEVICE inline float pseudo_huber(float delta, float r) {
    return pseudo_huber(static_cast<double>(delta), static_cast<double>(r));
}

} // namespace xsf
