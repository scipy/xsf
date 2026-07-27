#pragma once

#include "binom.h"
#include "config.h"
#include "hyp2f1.h"

namespace xsf {

namespace detail {

    template <typename T>
    XSF_HOST_DEVICE inline T eval_jacobi(double n, double alpha, double beta, T x) {
        double a, b, c, d;
        T g;

        if (alpha == -1 && std::abs(beta) == 1) {
            if (n == 0) {
                return 1.0;
            } else if (n == 1) {
                return 0.5 * (1.0 + beta) * (x - 1.0);
            } else if (n > 1) {
                return ((n + beta) / (2.0 * n)) * (x - 1.0) * eval_jacobi(n - 1.0, 1.0, beta, x);
            }
        }

        d = binom(n + alpha, n);
        a = -n;
        b = n + alpha + beta + 1.0;
        c = alpha + 1.0;
        g = 0.5 * (1.0 - x);
        return d * hyp2f1(a, b, c, g);
    }

    XSF_HOST_DEVICE inline double eval_jacobi_l(std::ptrdiff_t n, double alpha, double beta, double x) {
        std::ptrdiff_t kk;
        double p, d;
        double k, t;

        if (n < 0) {
            return eval_jacobi(n, alpha, beta, x);
        } else if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return 0.5 * (2.0 * (alpha + 1.0) + (alpha + beta + 2.0) * (x - 1.0));
        } else if (alpha == -1 && std::abs(beta) == 1) {
            return ((n + beta) / (2.0 * n)) * (x - 1.0) * eval_jacobi(n - 1.0, 1.0, beta, x);
        } else {
            d = (alpha + beta + 2.0) * (x - 1.0) / (2.0 * (alpha + 1.0));
            p = d + 1.0;
            for (kk = 0; kk < n - 1; kk++) {
                k = kk + 1.0;
                t = 2.0 * k + alpha + beta;
                d = ((t * (t + 1.0) * (t + 2.0)) * (x - 1.0) * p + 2.0 * k * (k + beta) * (t + 2.0) * d) /
                    (2.0 * (k + alpha + 1.0) * (k + alpha + beta + 1.0) * t);
                p = d + p;
            }
            return binom(n + alpha, n) * p;
        }
    }

} // namespace detail

// Jacobi

XSF_HOST_DEVICE inline double eval_jacobi(double n, double alpha, double beta, double x) {
    return detail::eval_jacobi(n, alpha, beta, x);
}

XSF_HOST_DEVICE inline float eval_jacobi(float n, float alpha, float beta, float x) {
    return detail::eval_jacobi(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(beta), static_cast<double>(x)
    );
}

XSF_HOST_DEVICE inline std::complex<double> eval_jacobi(double n, double alpha, double beta, std::complex<double> x) {
    return detail::eval_jacobi(n, alpha, beta, x);
}

XSF_HOST_DEVICE inline std::complex<float> eval_jacobi(float n, float alpha, float beta, std::complex<float> x) {
    return static_cast<std::complex<float>>(detail::eval_jacobi(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(beta),
        static_cast<std::complex<double>>(x)
    ));
}

XSF_HOST_DEVICE inline double eval_jacobi(std::ptrdiff_t n, double alpha, double beta, double x) {
    return detail::eval_jacobi_l(n, alpha, beta, x);
}

XSF_HOST_DEVICE inline float eval_jacobi(std::ptrdiff_t n, float alpha, float beta, float x) {
    return detail::eval_jacobi_l(n, static_cast<double>(alpha), static_cast<double>(beta), static_cast<double>(x));
}

// Shifted Jacobi

XSF_HOST_DEVICE inline double eval_sh_jacobi(double n, double p, double q, double x) {
    return detail::eval_jacobi(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

XSF_HOST_DEVICE inline double eval_sh_jacobi(std::ptrdiff_t n, double p, double q, double x) {
    return detail::eval_jacobi_l(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

XSF_HOST_DEVICE inline float eval_sh_jacobi(float n, float p, float q, float x) {
    return detail::eval_jacobi(
               static_cast<double>(n), static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
               2.0 * static_cast<double>(x) - 1.0
           ) /
           binom(2.0 * static_cast<double>(n) + static_cast<double>(p) - 1.0, static_cast<double>(n));
}

XSF_HOST_DEVICE inline float eval_sh_jacobi(std::ptrdiff_t n, float p, float q, float x) {
    return detail::eval_jacobi_l(
               n, static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
               2.0 * static_cast<double>(x) - 1.0
           ) /
           binom(2.0 * n + static_cast<double>(p) - 1.0, n);
}

XSF_HOST_DEVICE inline std::complex<double> eval_sh_jacobi(double n, double p, double q, std::complex<double> x) {
    return detail::eval_jacobi(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

XSF_HOST_DEVICE inline std::complex<float> eval_sh_jacobi(float n, float p, float q, std::complex<float> x) {
    return static_cast<std::complex<float>>(
        detail::eval_jacobi(
            static_cast<double>(n), static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
            2.0 * static_cast<std::complex<double>>(x) - 1.0
        ) /
        binom(2.0 * static_cast<double>(n) + static_cast<double>(p) - 1.0, static_cast<double>(n))
    );
}

} // namespace xsf
