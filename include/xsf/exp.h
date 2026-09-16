#pragma once

#include "config.h"
#include "trig.h"

#include "cephes/exp10.h"
#include "cephes/exp2.h"
#include "cephes/unity.h"

namespace xsf {

inline double expm1(double x) { return cephes::expm1(x); }

inline float expm1(float x) { return expm1(static_cast<double>(x)); }

// cexpm1(z) = cexp(z) - 1
//
// The imaginary part of this is easily computed via exp(z.real)*sin(z.imag)
// The real part is difficult to compute when there is cancellation e.g. when
// z.real = -log(cos(z.imag)).  There isn't a way around this problem  that
// doesn't involve computing exp(z.real) and/or cos(z.imag) to higher
// precision.
inline std::complex<double> expm1(std::complex<double> z) {
    if (!std::isfinite(std::real(z)) || !std::isfinite(std::imag(z))) {
        // libstdc++'s std::exp(complex) is implemented via std::polar and
        // does not follow C99 cexp semantics for non-finite inputs (e.g.
        // it returns (inf, nan) for exp(inf+0i) instead of (inf, 0), and
        // it aborts under _GLIBCXX_ASSERTIONS when a NaN reaches polar).
        // Compute the non-finite branches directly.
        double zr = std::real(z);
        double zi = std::imag(z);
        constexpr double inf = std::numeric_limits<double>::infinity();
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        if (std::isnan(zr)) {
            return std::complex<double>{nan, nan};
        }
        if (std::isinf(zr)) {
            if (zr < 0) {
                // exp(-inf + i*y) = 0 for any y (C99 conventions), so
                // expm1 = -1 + 0i regardless of zi.
                return std::complex<double>{-1.0, 0.0};
            }
            // zr == +inf
            if (zi == 0.0) {
                return std::complex<double>{inf, zi};
            }
            if (std::isfinite(zi)) {
                return std::complex<double>{inf * std::cos(zi), inf * std::sin(zi)};
            }
            // zi is +/-inf or NaN: magnitude is inf, phase undefined.
            return std::complex<double>{inf, nan};
        }
        // zr finite, zi non-finite (inf or NaN).
        return std::complex<double>{nan, nan};
    }

    double x;
    double ezr = 0;
    if (std::real(z) <= -40) {
        x = -1.0;
    } else {
        ezr = expm1(std::real(z));
        x = ezr * std::cos(std::imag(z)) + cosm1(std::imag(z));
    }

    // don't compute exp(zr) too, unless necessary
    double y;
    if (std::real(z) > -1.0) {
        y = (ezr + 1.0) * sin(std::imag(z));
    } else {
        y = exp(std::real(z)) * sin(std::imag(z));
    }

    return std::complex<double>{x, y};
}

inline std::complex<float> expm1(std::complex<float> z) {
    return static_cast<std::complex<float>>(expm1(static_cast<std::complex<double>>(z)));
}

inline double exp2(double x) { return cephes::exp2(x); }

inline float exp2(float x) { return exp2(static_cast<double>(x)); }

inline double exp10(double x) { return cephes::exp10(x); }

inline float exp10(float x) { return exp10(static_cast<double>(x)); }

} // namespace xsf
