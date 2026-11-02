#pragma once

#include "../orthogonal_eval.h"
#include "specfun.h"

namespace xsf {

namespace detail {

    template <typename T>
    inline T eval_genlaguerre(double n, double alpha, T x) {
        if (alpha <= -1) {
            set_error("eval_genlaguerre", SF_ERROR_DOMAIN, "polynomial defined only for alpha > -1");
            return std::numeric_limits<double>::quiet_NaN();
        }

        const double d = binom(n + alpha, n);
        return d * hyp1f1(-n, alpha + 1.0, x);
    }

} // namespace detail

// Generalized Laguerre

inline double eval_genlaguerre(double n, double alpha, double x) { return detail::eval_genlaguerre(n, alpha, x); }

inline float eval_genlaguerre(float n, float alpha, float x) {
    return detail::eval_genlaguerre(static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(x));
}

inline std::complex<double> eval_genlaguerre(double n, double alpha, std::complex<double> x) {
    return detail::eval_genlaguerre(n, alpha, x);
}

inline std::complex<float> eval_genlaguerre(float n, float alpha, std::complex<float> x) {
    return static_cast<std::complex<float>>(detail::eval_genlaguerre(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<std::complex<double>>(x)
    ));
}

// Laguerre

inline double eval_laguerre(double n, double x) { return detail::eval_genlaguerre(n, 0.0, x); }

inline float eval_laguerre(float n, float x) {
    return detail::eval_genlaguerre(static_cast<double>(n), 0.0, static_cast<double>(x));
}

inline std::complex<double> eval_laguerre(double n, std::complex<double> x) {
    return detail::eval_genlaguerre(n, 0.0, x);
}

inline std::complex<float> eval_laguerre(float n, std::complex<float> x) {
    return static_cast<std::complex<float>>(
        detail::eval_genlaguerre(static_cast<double>(n), 0.0, static_cast<std::complex<double>>(x))
    );
}

} // namespace xsf
