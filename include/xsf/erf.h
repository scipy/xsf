#pragma once

#include "cephes/ndtr.h"
#include "config.h"
#include "faddeeva.h"

namespace xsf {

XSF_HOST_DEVICE inline double erf(double x) { return cephes::erf(x); }

XSF_HOST_DEVICE inline float erf(float x) { return erf(static_cast<double>(x)); }

XSF_HOST_DEVICE inline cxx::complex<double> erf(cxx::complex<double> z) { return Faddeeva::erf(z); }

XSF_HOST_DEVICE inline cxx::complex<float> erf(cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(erf(static_cast<cxx::complex<double>>(x)));
}

XSF_HOST_DEVICE inline double erfc(double x) { return cephes::erfc(x); }

XSF_HOST_DEVICE inline float erfc(float x) { return erfc(static_cast<double>(x)); }

XSF_HOST_DEVICE inline cxx::complex<double> erfc(cxx::complex<double> z) { return Faddeeva::erfc(z); }

XSF_HOST_DEVICE inline cxx::complex<float> erfc(cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(erfc(static_cast<cxx::complex<double>>(x)));
}

XSF_HOST_DEVICE inline double erfcx(double x) { return Faddeeva::erfcx(x); }

XSF_HOST_DEVICE inline float erfcx(float x) { return erfcx(static_cast<double>(x)); }

XSF_HOST_DEVICE inline cxx::complex<double> erfcx(cxx::complex<double> z) { return Faddeeva::erfcx(z); }

XSF_HOST_DEVICE inline cxx::complex<float> erfcx(cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(erfcx(static_cast<cxx::complex<double>>(x)));
}

XSF_HOST_DEVICE inline double erfi(double x) { return Faddeeva::erfi(x); }

XSF_HOST_DEVICE inline float erfi(float x) { return erfi(static_cast<double>(x)); }

XSF_HOST_DEVICE inline cxx::complex<double> erfi(cxx::complex<double> z) { return Faddeeva::erfi(z); }

XSF_HOST_DEVICE inline cxx::complex<float> erfi(cxx::complex<float> z) {
    return static_cast<cxx::complex<float>>(erfi(static_cast<cxx::complex<double>>(z)));
}

XSF_HOST_DEVICE inline double voigt_profile(double x, double sigma, double gamma) {
    constexpr double INV_SQRT_2 = 0.707106781186547524401;
    constexpr double SQRT_2PI = 2.5066282746310002416123552393401042;

    if (sigma == 0) {
        if (gamma == 0) {
            if (cxx::isnan(x))
                return x;
            if (x == 0)
                return cxx::numeric_limits<double>::infinity();
            return 0;
        }
        return gamma / M_PI / (x * x + gamma * gamma);
    }
    if (gamma == 0) {
        return 1 / SQRT_2PI / sigma * exp(-(x / sigma) * (x / sigma) / 2);
    }

    double zreal = x / sigma * INV_SQRT_2;
    double zimag = gamma / sigma * INV_SQRT_2;
    cxx::complex<double> z(zreal, zimag);
    cxx::complex<double> w = Faddeeva::w(z);
    return w.real() / sigma / SQRT_2PI;
}

XSF_HOST_DEVICE inline float voigt_profile(float x, float sigma, float gamma) {
    return voigt_profile(static_cast<double>(x), static_cast<double>(sigma), static_cast<double>(gamma));
}

XSF_HOST_DEVICE inline cxx::complex<double> wofz(cxx::complex<double> z) { return Faddeeva::w(z); }

XSF_HOST_DEVICE inline cxx::complex<float> wofz(cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(wofz(static_cast<cxx::complex<double>>(x)));
}

XSF_HOST_DEVICE inline double dawsn(double x) { return Faddeeva::Dawson(x); }

XSF_HOST_DEVICE inline float dawsn(float x) { return dawsn(static_cast<double>(x)); }

XSF_HOST_DEVICE inline cxx::complex<double> dawsn(cxx::complex<double> z) { return Faddeeva::Dawson(z); }

XSF_HOST_DEVICE inline cxx::complex<float> dawsn(cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(dawsn(static_cast<cxx::complex<double>>(x)));
}

} // namespace xsf
