// Translated from Cython to C++ by the xsf developers in 2026.
// Original implementation by Warren Weckesser.
#pragma once

#include "config.h"

namespace xsf {

XSF_HOST_DEVICE inline double boxcox(double x, double lmbda) {
    /* If lmbda << 1 and log(x) < 1.0, the lmbda*log(x) product can lose
     * precision, furthermore, expm1(x) == x for x < eps.
     * For doubles, the range of log is -744.44 to +709.78, with eps being
     * the smallest value produced.  This range means that we will have
     * abs(lmbda)*log(x) < eps whenever abs(lmbda) <= eps/-log(min double)
     * which is ~2.98e-19.
     */
    if (std::abs(lmbda) < 1e-19) {
        return std::log(x);
    } else if (lmbda * std::log(x) < 709.78) {
        return std::expm1(lmbda * std::log(x)) / lmbda;
    } else {
        return std::copysign(1.0, lmbda) * std::exp(lmbda * std::log(x) - std::log(std::abs(lmbda))) - 1 / lmbda;
    }
}

XSF_HOST_DEVICE inline float boxcox(float x, float lmbda) {
    return (boxcox(static_cast<double>(x), static_cast<double>(lmbda)));
}

XSF_HOST_DEVICE inline double boxcox1p(double x, double lmbda) {
    /* The argument given above in boxcox applies here with the modification
     * that the smallest value produced by log1p is the minimum representable
     * value, rather than eps.  The second condition here prevents underflow
     * when log1p(x) is < eps.
     */
    double lgx = std::log1p(x);
    if (std::abs(lmbda) < 1e-19 || (std::abs(lgx) < 1e-289 && std::abs(lmbda) < 1e273)) {
        return lgx;
    } else if (lmbda * lgx < 709.78) {
        return std::expm1(lmbda * lgx) / lmbda;
    } else {
        return std::copysign(1.0, lmbda) * std::exp(lmbda * lgx - std::log(std::abs(lmbda))) - 1 / lmbda;
    }
}

XSF_HOST_DEVICE inline float boxcox1p(float x, float lmbda) {
    return (boxcox1p(static_cast<double>(x), static_cast<double>(lmbda)));
}

XSF_HOST_DEVICE inline double inv_boxcox(double x, double lmbda) {
    if (lmbda == 0) {
        return std::exp(x);
    } else if (lmbda * x < 1.79e308) {
        return std::exp(std::log1p(lmbda * x) / lmbda);
    } else {
        return std::exp((std::log(std::copysign(1.0, lmbda) * (x + 1 / lmbda)) + std::log(std::abs(lmbda))) / lmbda);
    }
}

XSF_HOST_DEVICE inline float inv_boxcox(float x, float lmbda) {
    return (inv_boxcox(static_cast<double>(x), static_cast<double>(lmbda)));
}

XSF_HOST_DEVICE inline double inv_boxcox1p(double x, double lmbda) {
    if (lmbda == 0) {
        return std::expm1(x);
    } else if (std::abs(lmbda * x) < 1e-154) {
        return x;
    } else if (lmbda * x < 1.79e308) {
        return std::expm1(std::log1p(lmbda * x) / lmbda);
    } else {
        return std::expm1((std::log(std::copysign(1.0, lmbda) * (x + 1 / lmbda)) + std::log(std::abs(lmbda))) / lmbda);
    }
}

XSF_HOST_DEVICE inline float inv_boxcox1p(float x, float lmbda) {
    return (inv_boxcox1p(static_cast<double>(x), static_cast<double>(lmbda)));
}

} // namespace xsf
