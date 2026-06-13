#pragma once

#include "cephes/gamma.h"
#include "cephes/igam.h"
#include "cephes/igami.h"
#include "loggamma.h"

namespace xsf {

template <typename T>
XSF_HOST_DEVICE T gamma(T x) {
    return cephes::Gamma(x);
}

XSF_HOST_DEVICE inline double gammainc(double a, double x) { return cephes::igam(a, x); }

XSF_HOST_DEVICE inline float gammainc(float a, float x) {
    return gammainc(static_cast<double>(a), static_cast<double>(x));
}

XSF_HOST_DEVICE inline double gammaincinv(double a, double p) { return cephes::igami(a, p); }

XSF_HOST_DEVICE inline float gammaincinv(float a, float p) {
    return gammaincinv(static_cast<double>(a), static_cast<double>(p));
}

XSF_HOST_DEVICE inline double gammaincc(double a, double x) { return cephes::igamc(a, x); }

XSF_HOST_DEVICE inline float gammaincc(float a, float x) {
    return gammaincc(static_cast<double>(a), static_cast<double>(x));
}

XSF_HOST_DEVICE inline double gammainccinv(double a, double p) { return cephes::igamci(a, p); }

XSF_HOST_DEVICE inline float gammainccinv(float a, float p) {
    return gammainccinv(static_cast<double>(a), static_cast<double>(p));
}

XSF_HOST_DEVICE inline double gammaln(double x) { return cephes::lgam(x); }

XSF_HOST_DEVICE inline float gammaln(float x) { return gammaln(static_cast<double>(x)); }

XSF_HOST_DEVICE inline double gammasgn(double x) { return cephes::gammasgn(x); }

XSF_HOST_DEVICE inline float gammasgn(float x) { return gammasgn(static_cast<double>(x)); }

XSF_HOST_DEVICE inline std::complex<double> gamma(std::complex<double> z) {
    // Guard against NaN/Inf inputs: std::exp(complex) is implemented in
    // libstdc++ as std::polar(std::exp(re), im), and std::polar asserts
    // __rho >= 0 under _GLIBCXX_ASSERTIONS -- which is false for NaN.
    if (!std::isfinite(z.real()) || !std::isfinite(z.imag())) {
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    // Compute Gamma(z) using loggamma.
    if (z.real() <= 0 && z == std::floor(z.real())) {
        // Poles
        set_error("gamma", SF_ERROR_SINGULAR, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }

    if (z.real() <= -std::ldexp(1.0, std::numeric_limits<double>::digits)) {
        return {0.0, std::copysign(0.0, z.imag())};
    }

    std::complex<double> lg = loggamma(z);
    if (lg.real() == -std::numeric_limits<double>::infinity()) {
        return {0.0, std::copysign(0.0, z.imag())};
    }
    const double max = std::numeric_limits<double>::max();
    if (lg.real() > std::log(max) && z.imag() == 0.0) {
        return {std::numeric_limits<double>::infinity(), std::copysign(0.0, z.imag())};
    }
    if (lg.real() > std::log(max) && z.real() > std::sqrt(max) && std::abs(z.imag()) > std::sqrt(max)) {
        return {
            -std::numeric_limits<double>::infinity(), std::copysign(std::numeric_limits<double>::infinity(), z.imag())
        };
    }
    return std::exp(lg);
}

XSF_HOST_DEVICE inline std::complex<float> gamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(gamma(static_cast<std::complex<double>>(z)));
}

template <typename T>
T gamma_ratio(T a, T b) {
    return std::tgamma(a) / std::tgamma(b);
}

} // namespace xsf
