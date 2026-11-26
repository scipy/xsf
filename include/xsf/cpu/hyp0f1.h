// Translated from Cython to C++ by the xsf developers in 2026.
#pragma once

#include "../config.h"
#include "../gamma.h"
#include "../log.h"
#include "../trig.h"
#include "bessel.h"

namespace xsf {
namespace detail {

    template <typename T>
    T hyp0f1_asy(double v, T z) {
        /* Asymptotic expansion for I_{v-1}(2*sqrt(z)) * Gamma(v)
         * for v -> +infinity, with z not near (-inf, -(v-1)**2/4].
         *
         * Based off DLMF 10.41
         */
        T arg = std::sqrt(z);
        T v1 = std::abs(v - 1.0);
        T x = 2.0 * arg / v1;
        T p1 = std::sqrt(1.0 + x * x);
        T eta = p1 + std::log(x) - xsf::log1p(p1);
        T pp, p2, p4, p6, u1, u2, u3, u_corr_i, u_corr_k;
        T result;
        double gs;

        T arg_exp_i = -0.5 * std::log(p1);
        arg_exp_i -= 0.5 * std::log(2.0 * M_PI * v1);
        arg_exp_i += xsf::gammaln(v);
        gs = xsf::gammasgn(v);

        T arg_exp_k = arg_exp_i;
        arg_exp_i += v1 * eta;
        arg_exp_k -= v1 * eta;

        // Large-v asymptotic correction, DLMF 10.41.10.
        pp = 1.0 / p1;
        p2 = pp * pp;
        p4 = p2 * p2;
        p6 = p4 * p2;
        u1 = (3.0 - 5.0 * p2) * pp / 24.0;
        u2 = (81.0 - 462.0 * p2 + 385.0 * p4) * p2 / 1152.0;
        u3 = (30375.0 - 369603.0 * p2 + 765765.0 * p4 - 425425.0 * p6) * pp * p2 / 414720.0;
        u_corr_i = 1.0 + u1 / v1 + u2 / (v1 * v1) + u3 / (v1 * v1 * v1);

        // DLMF 10.39.9: 0F1(;v;z) = Gamma(v) arg^(1-v) I_{v-1}(2 arg)
        result = std::exp(arg_exp_i + xsf::xlogy(T(1.0 - v), arg)) * gs * u_corr_i;
        if (v - 1.0 < 0.0) {
            // DLMF 10.27.2: I_{-v} = I_{v} + (2/pi) sin(pi*v) K_v.
            u_corr_k = 1.0 - u1 / v1 + u2 / (v1 * v1) - u3 / (v1 * v1 * v1);
            result += std::exp(arg_exp_k + xsf::xlogy(v1, arg)) * gs * 2.0 * xsf::sinpi(v1) * u_corr_k;
        }

        return result;
    }

} // namespace detail

inline double hyp0f1(double v, double z) {
    double arg, arg_exp, bess_val;

    if (std::isnan(v) || std::isnan(z)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Handle poles and zeros.
    if (v <= 0.0 && v == std::floor(v)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (z == 0.0 && v != 0.0) {
        return 1.0;
    }

    // Both v and z small: truncate the Taylor series at O(z**2).
    if (std::abs(z) < 1e-6 * (1.0 + std::abs(v))) {
        return 1.0 + z / v + z * z / (2.0 * v * (v + 1.0));
    }

    if (z > 0.0) {
        arg = std::sqrt(z);
        arg_exp = xsf::xlogy(1.0 - v, arg) + xsf::gammaln(v);
        bess_val = xsf::cyl_bessel_i(v - 1.0, 2.0 * arg);

        if (arg_exp > std::log(std::numeric_limits<double>::max()) || bess_val == 0.0 ||      // overflow
            arg_exp < std::log(std::numeric_limits<double>::min()) || std::isinf(bess_val)) { // underflow
            return detail::hyp0f1_asy(v, z);
        }
        return std::exp(arg_exp) * xsf::gammasgn(v) * bess_val;
    }

    arg = std::sqrt(-z);
    double arg_pow = std::pow(arg, 1.0 - v);
    double gamma_v = xsf::gamma(v);
    bess_val = xsf::cyl_bessel_j(v - 1.0, 2.0 * arg);
    double prefactor = arg_pow * gamma_v;
    if (std::isnormal(arg_pow) && std::isnormal(gamma_v) && std::isnormal(prefactor) && std::isnormal(bess_val)) {
        return prefactor * bess_val;
    }

    // The power, Gamma(v) or their product over- or underflows.
    if (std::isnormal(bess_val)) {
        arg_exp = xsf::xlogy(1.0 - v, arg) + xsf::gammaln(v) + std::log(std::abs(bess_val));
        return std::copysign(std::exp(arg_exp), bess_val) * xsf::gammasgn(v);
    }

    // The Bessel function under- (v > 1) or overflows (v < 1) below its turning point: Debye's expansion
    // (DLMF 10.19.3), i.e. hyp0f1_asy at z < 0.
    if ((v > 1.0 || std::isinf(bess_val)) && 2.0 * arg < std::abs(v - 1.0)) {
        return detail::hyp0f1_asy(v, std::complex<double>(z)).real();
    }
    return prefactor * bess_val;
}

inline float hyp0f1(float v, float z) { return hyp0f1(static_cast<double>(v), static_cast<double>(z)); }

inline std::complex<double> hyp0f1(double v, std::complex<double> z) {
    std::complex<double> arg, s, t1, t2;
    std::complex<double> r;

    if (std::isnan(v) || std::isnan(z.real()) || std::isnan(z.imag())) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan, nan};
    }

    // Handle poles and zeros.
    if (v <= 0.0 && v == std::floor(v)) {
        return {std::numeric_limits<double>::quiet_NaN(), 0.0};
    }
    if (z.real() == 0.0 && z.imag() == 0.0 && v != 0.0) {
        return {1.0, 0.0};
    }

    // Both v and z small: truncate the Taylor series at O(z**2).
    if (std::abs(z) < 1e-6 * (1.0 + std::abs(v))) {
        // Need to do computations in this order, for otherwise $v\approx -z \ll 1$.
        // It can lose precision (as was reported for 32-bit linux, see gh-6365)
        t1 = 1.0 + z / v;
        t2 = z * z / (2.0 * v * (v + 1.0));
        return t1 + t2;
    }

    if (z.real() > 0.0) {
        arg = std::sqrt(z);
        s = 2.0 * arg;
        r = xsf::cyl_bessel_i(v - 1.0, s);
    } else {
        arg = std::sqrt(-z);
        s = 2.0 * arg;
        r = xsf::cyl_bessel_j(v - 1.0, s);
    }

    t1 = r * xsf::gamma(v) * std::pow(arg, 1.0 - v);
    if (t1 != 0.0 && std::isfinite(t1.real()) && std::isfinite(t1.imag())) {
        return t1;
    }

    // The Bessel function over- or underflows.
    if (v > 1.0 && (r == 0.0 || !std::isfinite(std::abs(r)))) {
        return detail::hyp0f1_asy(v, z);
    }

    // Gamma(v) or the power over- or underflows.
    return std::exp(xsf::xlogy(std::complex<double>(1.0 - v), arg) + xsf::gammaln(v) + std::log(r)) * xsf::gammasgn(v);
}

inline std::complex<float> hyp0f1(float v, std::complex<float> z) {
    return static_cast<std::complex<float>>(hyp0f1(static_cast<double>(v), static_cast<std::complex<double>>(z)));
}

} // namespace xsf
