#pragma once

#include "binom.h"
#include "config.h"
#include "hyp2f1.h"

namespace xsf {

// Jacobi

template <typename T>
XSF_HOST_DEVICE inline T eval_jacobi(double n, double alpha, double beta, T x) {
    double a, b, c, d;
    T g;

    if (alpha == -1 and std::abs(beta) == 1) {
        if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return 0.5 * (1.0 + beta) * (x - 1.0);
        } else if (n > 1) {
            return ((n + beta) / (2.0 * n)) * (x - 1) * eval_jacobi(n - 1, 1, beta, x);
        }
    }

    d = binom(n + alpha, n);
    a = -n;
    b = n + alpha + beta + 1;
    c = alpha + 1;
    g = 0.5 * (1 - x);
    return d * hyp2f1(a, b, c, g);
}

XSF_HOST_DEVICE inline double eval_jacobi_l(int n, double alpha, double beta, double x) {
    int kk;
    double p, d;
    double k, t;

    if (n < 0) {
        return eval_jacobi(n, alpha, beta, x);
    } else if (n == 0) {
        return 1.0;
    } else if (n == 1) {
        return 0.5 * (2 * (alpha + 1) + (alpha + beta + 2) * (x - 1));
    } else if (alpha == -1 and std::abs(beta) == 1) {
        return ((n + beta) / (2.0 * n)) * (x - 1) * eval_jacobi(n - 1, 1, beta, x);
    } else {
        d = (alpha + beta + 2) * (x - 1) / (2 * (alpha + 1));
        p = d + 1;
        for (kk = 0; kk < n - 1; kk++) {
            k = kk + 1.0;
            t = 2 * k + alpha + beta;
            d = ((t * (t + 1) * (t + 2)) * (x - 1) * p + 2 * k * (k + beta) * (t + 2) * d) /
                (2 * (k + alpha + 1) * (k + alpha + beta + 1) * t);
            p = d + p;
        }
        return binom(n + alpha, n) * p;
    }
}

template <typename T>
XSF_HOST_DEVICE inline float eval_jacobi(float n, float alpha, float beta, T x) {
    return eval_jacobi(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(beta), static_cast<double>(x)
    );
}

XSF_HOST_DEVICE inline float eval_jacobi_l(int n, float alpha, float beta, float x) {
    return eval_jacobi_l(n, static_cast<double>(alpha), static_cast<double>(beta), static_cast<double>(x));
}

// Shifted Jacobi

template <typename T>
XSF_HOST_DEVICE inline T eval_sh_jacobi(double n, double p, double q, T x) {
    return eval_jacobi(n, p - q, q - 1, 2 * x - 1) / binom(2 * n + p - 1, n);
}

XSF_HOST_DEVICE inline double eval_sh_jacobi_l(int n, double p, double q, double x) {
    return eval_jacobi_l(n, p - q, q - 1, 2 * x - 1) / binom(2 * n + p - 1, n);
}

template <typename T>
XSF_HOST_DEVICE inline T eval_sh_jacobi(float n, float p, float q, T x) {
    return eval_sh_jacobi(
        static_cast<double>(n), static_cast<double>(p), static_cast<double>(q), static_cast<double>(x)
    );
}

XSF_HOST_DEVICE inline float eval_sh_jacobi_l(int n, float p, float q, float x) {
    return eval_sh_jacobi_l(n, static_cast<double>(p), static_cast<double>(q), static_cast<double>(x));
}

} // namespace xsf
