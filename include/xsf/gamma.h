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

XSF_HOST_DEVICE inline cxx::complex<double> gamma(cxx::complex<double> z) {
    // Guard against NaN/Inf inputs: std::exp(complex) is implemented in
    // libstdc++ as std::polar(std::exp(re), im), and std::polar asserts
    // __rho >= 0 under _GLIBCXX_ASSERTIONS -- which is false for NaN.
    if (!cxx::isfinite(z.real()) || !cxx::isfinite(z.imag())) {
        return {cxx::numeric_limits<double>::quiet_NaN(), cxx::numeric_limits<double>::quiet_NaN()};
    }
    // Compute Gamma(z) using loggamma.
    if (z.real() <= 0 && z == cxx::floor(z.real())) {
        // Gamma poles at non-positive integers.
        set_error("gamma", SF_ERROR_SINGULAR, NULL);
        return {cxx::numeric_limits<double>::quiet_NaN(), cxx::numeric_limits<double>::quiet_NaN()};
    }

    if (z.real() <= -cxx::ldexp(1.0, cxx::numeric_limits<double>::digits)) {
        // For real(z) <= -2**53, every representable real part has even-integer
        // parity, and Gamma(z) underflows to signed zero.
        return {0.0, cxx::copysign(0.0, z.imag())};
    }

    cxx::complex<double> lg = loggamma(z);
    if (lg.real() == -cxx::numeric_limits<double>::infinity()) {
        return {0.0, cxx::copysign(0.0, z.imag())};
    }
    const double max = cxx::numeric_limits<double>::max();
    if (lg.real() > cxx::log(max) && z.imag() == 0.0) {
        return {cxx::numeric_limits<double>::infinity(), cxx::copysign(0.0, z.imag())};
    }
    if (lg.real() > cxx::log(max) && z.real() > cxx::sqrt(max) && cxx::abs(z.imag()) > cxx::sqrt(max)) {
        // Avoid std::exp(complex) overflow; the quadrant follows sign(imag(z)).
        return {
            -cxx::numeric_limits<double>::infinity(), cxx::copysign(cxx::numeric_limits<double>::infinity(), z.imag())
        };
    }
    return cxx::exp(lg);
}

XSF_HOST_DEVICE inline cxx::complex<float> gamma(cxx::complex<float> z) {
    return static_cast<cxx::complex<float>>(gamma(static_cast<cxx::complex<double>>(z)));
}

template <typename T>
T gamma_ratio(T a, T b) {
    return std::tgamma(a) / std::tgamma(b);
}

} // namespace xsf
