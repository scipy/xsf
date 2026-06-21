// Translated from Cython to C++ by the xsf developers in 2026.
#pragma once

#include "bessel.h"
#include "config.h"
#include "gamma.h"
#include "log.h"
#include "trig.h"

namespace xsf {
namespace detail {

    inline double hyp0f1_asy(double v, double z) {
        /* Asymptotic expansion for I_{v-1}(2*sqrt(z)) * Gamma(v)
         * for real z > 0 and v -> +infinity.
         *
         * Based off DLMF 10.41
         */
        double arg = std::sqrt(z);
        double v1 = std::abs(v - 1.0);
        double x = 2.0 * arg / v1;
        double p1 = std::sqrt(1.0 + x * x);
        double eta = p1 + std::log(x) - std::log1p(p1);
        double pp, p2, p4, p6, u1, u2, u3, u_corr_i, u_corr_k;
        double result, gs;

        double arg_exp_i = -0.5 * std::log(p1);
        arg_exp_i -= 0.5 * std::log(2.0 * M_PI * v1);
        arg_exp_i += xsf::gammaln(v);
        gs = xsf::gammasgn(v);

        double arg_exp_k = arg_exp_i;
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

        result = std::exp(arg_exp_i - xsf::xlogy(v1, arg)) * gs * u_corr_i;
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
    return std::pow(arg, 1.0 - v) * xsf::gamma(v) * xsf::cyl_bessel_j(v - 1.0, 2.0 * arg);
}

inline float hyp0f1(float v, float z) { return hyp0f1(static_cast<double>(v), static_cast<double>(z)); }

inline std::complex<double> hyp0f1(double v, std::complex<double> z) {
    std::complex<double> arg, s, t1, t2;
    std::complex<double> r;

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

    return r * xsf::gamma(v) * std::pow(arg, 1.0 - v);
}

inline std::complex<float> hyp0f1(float v, std::complex<float> z) {
    return static_cast<std::complex<float>>(hyp0f1(static_cast<double>(v), static_cast<std::complex<double>>(z)));
}

} // namespace xsf
