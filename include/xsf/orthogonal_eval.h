#pragma once

#include "binom.h"
#include "config.h"
#include "error.h"
#include "hyp2f1.h"

namespace xsf {

namespace detail {

    template <typename T>
    XSF_HOST_DEVICE inline T eval_jacobi(double n, double alpha, double beta, T x) {
        double a, b, c, d;
        T g;

        if (alpha == -1 && cxx::abs(beta) == 1) {
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

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_jacobi_l(Int n, double alpha, double beta, double x) {
        Int kk;
        double p, d;
        double k, t;

        if (n < 0) {
            return eval_jacobi(n, alpha, beta, x);
        } else if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return 0.5 * (2.0 * (alpha + 1.0) + (alpha + beta + 2.0) * (x - 1.0));
        } else if (alpha == -1 && cxx::abs(beta) == 1) {
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

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_hermitenorm(Int n, double x) {
        if (cxx::isnan(x)) {
            return x;
        }

        if (n < 0) {
            set_error("eval_hermitenorm", SF_ERROR_DOMAIN, "polynomial only defined for nonnegative n");
            return cxx::numeric_limits<double>::quiet_NaN();
        } else if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return x;
        }

        double y3 = 0.0;
        double y2 = 1.0;
        for (Int k = n; k > 1; --k) {
            const double y1 = x * y2 - k * y3;
            y3 = y2;
            y2 = y1;
        }
        return x * y2 - y3;
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_hermite(Int n, double x) {
        if (n < 0) {
            set_error("eval_hermite", SF_ERROR_DOMAIN, "polynomial only defined for nonnegative n");
            return cxx::numeric_limits<double>::quiet_NaN();
        }
        return eval_hermitenorm(n, cxx::sqrt(2.0) * x) * cxx::pow(2.0, n / 2.0);
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

XSF_HOST_DEVICE inline cxx::complex<double> eval_jacobi(double n, double alpha, double beta, cxx::complex<double> x) {
    return detail::eval_jacobi(n, alpha, beta, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_jacobi(float n, float alpha, float beta, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(detail::eval_jacobi(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(beta),
        static_cast<cxx::complex<double>>(x)
    ));
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_jacobi(Int n, double alpha, double beta, double x) {
    return detail::eval_jacobi_l(n, alpha, beta, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_jacobi(Int n, float alpha, float beta, float x) {
    return detail::eval_jacobi_l(n, static_cast<double>(alpha), static_cast<double>(beta), static_cast<double>(x));
}

// Shifted Jacobi

XSF_HOST_DEVICE inline double eval_sh_jacobi(double n, double p, double q, double x) {
    return detail::eval_jacobi(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_sh_jacobi(Int n, double p, double q, double x) {
    return detail::eval_jacobi_l(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

XSF_HOST_DEVICE inline float eval_sh_jacobi(float n, float p, float q, float x) {
    return detail::eval_jacobi(
               static_cast<double>(n), static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
               2.0 * static_cast<double>(x) - 1.0
           ) /
           binom(2.0 * static_cast<double>(n) + static_cast<double>(p) - 1.0, static_cast<double>(n));
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_sh_jacobi(Int n, float p, float q, float x) {
    return detail::eval_jacobi_l(
               n, static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
               2.0 * static_cast<double>(x) - 1.0
           ) /
           binom(2.0 * n + static_cast<double>(p) - 1.0, n);
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_sh_jacobi(double n, double p, double q, cxx::complex<double> x) {
    return detail::eval_jacobi(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_sh_jacobi(float n, float p, float q, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_jacobi(
            static_cast<double>(n), static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
            2.0 * static_cast<cxx::complex<double>>(x) - 1.0
        ) /
        binom(2.0 * static_cast<double>(n) + static_cast<double>(p) - 1.0, static_cast<double>(n))
    );
}

// Hermite (probabilist's)

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_hermitenorm(Int n, double x) {
    return detail::eval_hermitenorm(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_hermitenorm(Int n, float x) {
    return detail::eval_hermitenorm(n, static_cast<double>(x));
}

// Hermite (physicist's)

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_hermite(Int n, double x) {
    return detail::eval_hermite(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_hermite(Int n, float x) {
    return detail::eval_hermite(n, static_cast<double>(x));
}

} // namespace xsf
