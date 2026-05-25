#include "bessel.h"
#include "cephes/kolmogorov.h"
#include "config.h"
#include "gamma.h"

namespace xsf {
namespace cpu {

    inline double kolmogorov(double x) { return cephes::kolmogorov(x); }

    inline double kolmogc(double x) { return cephes::kolmogc(x); }

    inline double kolmogi(double x) { return cephes::kolmogi(x); }

    inline double kolmogci(double x) { return cephes::kolmogci(x); }

    inline double kolmogp(double x) { return cephes::kolmogp(x); }

    inline double smirnov(int n, double x) { return cephes::smirnov(n, x); }

    inline double smirnovc(int n, double x) { return cephes::smirnovc(n, x); }

    inline double smirnovi(int n, double x) { return cephes::smirnovi(n, x); }

    inline double smirnovci(int n, double x) { return cephes::smirnovci(n, x); }

    inline double smirnovp(int n, double x) { return cephes::smirnovp(n, x); }

    inline double cdf_cvm_inf(double x) {
        // CDF of the Cramer-von Mises test statistic (infinite sample limit).
        // Accurate for practical hypothesis testing but not expected to be accurate
        // for large values of x, e.g. x > 4, when the cdf is very close to 1.
        if (std::isnan(x)) {
            return x;
        }
        if (x <= 0) {
            return 0.0;
        }

        double total = 0.0;
        int k = 0;

        while (true) {
            double u = std::exp(gammaln(k + 0.5) - gammaln(k + 1.0)) / (std::pow(M_PI, 1.5) * std::sqrt(x));
            double y = 4.0 * k + 1.0;
            double q = (y * y) / (16.0 * x);
            double b = cyl_bessel_k(0.25, q);
            double z = u * std::sqrt(y) * std::exp(-q) * b;

            total += z;
            if (std::abs(z) < 1e-7) {
                break;
            }
            ++k;
        }

        return total;
    }

    // Helpers for psi1 computations
    XSF_HOST_DEVICE inline double psi1_G(double y) {
        // Helper implementing the G-type Bessel combination appearing in equation 1.10
        // in Csörgő, S. and Faraway, J. (1996).
        double z = y * y / 4.0;
        double b = cyl_bessel_k(0.25, z) + cyl_bessel_k(0.75, z);
        return std::exp(-z) * std::pow(y / 2.0, 1.5) * b / std::sqrt(M_PI);
    }

    XSF_HOST_DEVICE inline double psi1_H(double y) {
        // Helper implementing the H-type Bessel combination appearing in equation 1.10
        // in Csörgő, S. and Faraway, J. (1996).
        double z = y * y / 4.0;
        double b = 2.0 * cyl_bessel_k(0.25, z) + 3.0 * cyl_bessel_k(0.75, z) - cyl_bessel_k(1.25, z);
        return std::exp(-z) * std::pow(y / 2.0, 2.5) * b / std::sqrt(M_PI);
    }

    XSF_HOST_DEVICE inline double psi1_term(int k, double x) {
        // Compute the k-th term of the psi1 series expansion in equation 1.10
        // of Csörgő, S. and Faraway, J. (1996).
        double m = 2.0 * k + 1.0;
        double sx = 2.0 * std::sqrt(x);
        double y1 = std::pow(x, 0.75);
        double y2 = std::pow(x, 1.25);
        double gamma_kp1_2 = gamma(k + 0.5);
        double gamma_kp3_2 = gamma(k + 1.5);
        double e1 = m * gamma_kp1_2 * psi1_G((4.0 * k + 3.0) / sx) / (9.0 * y1);
        double e2 = gamma_kp1_2 * psi1_H((4.0 * k + 1.0) / sx) / (72.0 * y2);
        double e3 = 2.0 * (m + 2.0) * gamma_kp3_2 * psi1_H((4.0 * k + 5.0) / sx) / (12.0 * y2);
        double e4 = 7.0 * m * gamma_kp1_2 * psi1_G((4.0 * k + 1.0) / sx) / (144.0 * y1);
        double e5 = 7.0 * m * gamma_kp1_2 * psi1_G((4.0 * k + 5.0) / sx) / (144.0 * y1);
        return e1 + e2 + e3 + e4 + e5;
    }

    XSF_HOST_DEVICE inline double psi1_mod(double x) {
        // psi1 is defined in equation 1.10 in Csörgő, S. and Faraway, J. (1996).
        // This implements a modified version by excluding the term V(x) / 12
        // (here: cdf_cvm_inf(x) / 12) to avoid evaluating cdf_cvm_inf(x) twice in cdf_cvm.
        //
        // Implementation based on MAPLE code of Julian Faraway and R code of the
        // function pCvM in the package goftest (v1.1.1), permission granted
        // by Adrian Baddeley. Main difference in the implementation: the code
        // here keeps adding terms of the series until the terms are small enough.
        int k = 0;
        double tot = 0.0;
        while (true) {
            double gamma_kp1 = gamma(k + 1.0);
            double z = -psi1_term(k, x) / (M_PI * gamma_kp1);
            tot += z;
            if (std::abs(z) < 1e-7) {
                break;
            }
            ++k;
        }
        return tot;
    }

    XSF_HOST_DEVICE inline double cdf_cvm(double x, int n = -1) {
        // Calculate the CDF of the Cramér-von Mises statistic for a finite sample
        // size n. If n=-1, use the asymptotic CDF (n=inf).
        //
        // See equation 1.8 in Csörgő, S. and Faraway, J. (1996) for finite samples,
        // 1.2 for the asymptotic CDF.
        //
        // For finite n, the approximation is less accurate near the support boundaries
        // [1/(12*n), n/3] and for larger x values where the CDF is close to 1.
        // The returned value is clipped to [0, 1].
        if (n == -1) {
            return cdf_cvm_inf(x);
        }
        double lower = 1.0 / (12.0 * n);
        double upper = static_cast<double>(n) / 3.0;
        if (x <= lower) {
            return 0.0;
        }
        if (x >= upper) {
            return 1.0;
        }
        double cdf = cdf_cvm_inf(x) * (1.0 + lower) + psi1_mod(x) / n;
        return std::min(std::max(cdf, 0.0), 1.0);
    }

} // namespace cpu
} // namespace xsf
