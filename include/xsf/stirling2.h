/* Stirling numbers of the second kind
 *
 * SYNOPSIS: Stirling numbers of the second kind count the number of ways to
 *  make a partition of n distinct elements into k non-empty subsets.
 *
 * DESCRIPTION: n is the number of distinct elements of the set to be
 *  partitioned and k is the number of non-empty subsets. Values for n < 0 or
 *  k < 0 are interpreted as 0.
 *
 * NOTE: Original implementation by Lucas Roberts. Ported to xsf by SciPy
 *  developers.
 *
 * References:
 *  - N. M. Temme, "Asymptotic estimates of Stirling numbers", Studies in
 *    Applied Mathematics, 89(3):233-243, 1993.
 */

#pragma once

#include "binom.h"
#include "config.h"
#include "error.h"
#include "lambertw.h"

namespace xsf {
namespace detail {

    // Dynamic programming method for small n.
    // Only valid for n <= 50; the stack array is sized accordingly (max 25 elements).
    XSF_HOST_DEVICE inline double stirling2_dp(double n, double k) {
        if ((n == 0 && k == 0) || (n == 1 && k == 1)) {
            return 1.0;
        }
        if (k <= 0 || k > n || n < 0) {
            return 0.0;
        }

        int ni = static_cast<int>(n);
        int ki = static_cast<int>(k);
        int arraySize = ki <= ni - ki + 1 ? ki : ni - ki + 1;

        // Maximum arraySize for n <= 50 is 25. Guard against misuse.
        constexpr int max_array_size = 26;
        if (arraySize > max_array_size) {
            set_error("stirling2", SF_ERROR_OTHER, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }

        double curr[max_array_size];
        for (int i = 0; i < arraySize; i++) {
            curr[i] = 1.0;
        }

        if (ki <= ni - ki + 1) {
            for (int i = 1; i < ni - ki + 1; i++) {
                for (int j = 1; j < ki; j++) {
                    curr[j] = (j + 1) * curr[j] + curr[j - 1];
                    if (std::isinf(curr[j])) {
                        set_error("stirling2", SF_ERROR_OVERFLOW, NULL);
                        return std::numeric_limits<double>::infinity();
                    }
                }
            }
        } else {
            for (int i = 1; i < ki; i++) {
                for (int j = 1; j < ni - ki + 1; j++) {
                    curr[j] = (i + 1) * curr[j - 1] + curr[j];
                    if (std::isinf(curr[j])) {
                        set_error("stirling2", SF_ERROR_OVERFLOW, NULL);
                        return std::numeric_limits<double>::infinity();
                    }
                }
            }
        }
        return curr[arraySize - 1];
    }

    // Second-order Temme asymptotic approximation for large n.
    // Reference: Temme (1993).
    XSF_HOST_DEVICE inline double stirling2_temme(double n, double k) {
        if ((n == k && n >= 0) || (n > 0 && k == 1)) {
            return 1.0;
        }
        if (k <= 0 || k > n || n < 0) {
            return 0.0;
        }

        double mu = k / n;
        double d = std::exp(-1.0 / mu) / mu;
        std::complex<double> delta(-d, 0.0);

        // lambertw returns complex; we only need the real part (branch k=0)
        std::complex<double> lwv = xsf::lambertw(delta, 0, 1e-8);
        double x0 = lwv.real() + 1.0 / mu;
        double t0 = 1.0 / mu - 1.0;
        double F = std::sqrt(t0 / ((1.0 + t0) * (x0 - t0)));
        double A = -n * std::log(x0) + k * std::log(std::exp(x0) - 1.0) - k * t0 + (n - k) * std::log(t0);

        // F1 correction term (Horner scheme applied to numerator)
        double xt = x0 * t0;
        double t0power3 = t0 * t0 * t0;
        double num = -2.0 * x0 * x0 * x0;
        num += ((t0 + 2.0) * t0 + 2.0) * (2.0 * t0power3);
        num += (-6.0 * t0power3 + (8.0 * t0 - 6.0 * x0 - 5.0) * xt + ((2.0 * x0 + 1.0) * x0 + 3.0) * x0) * xt;
        double denom = 24.0 * F * (1.0 + t0) * (1.0 + t0) * (x0 - t0) * (x0 - t0) * (x0 - t0) * (x0 - t0);
        double F1 = num / denom;

        return std::exp(A) * std::pow(k, n - k) * xsf::binom(n, k) * (F - F1 / k);
    }

    XSF_HOST_DEVICE inline double stirling2(double n, double k) {
        if (n <= 50) {
            return detail::stirling2_dp(n, k);
        }
        return detail::stirling2_temme(n, k);
    }
} // namespace detail

XSF_HOST_DEVICE inline double stirling2(double n, double k) { return detail::stirling2(n, k); }

XSF_HOST_DEVICE inline float stirling2(float n, float k) {
    return static_cast<float>(detail::stirling2(static_cast<double>(n), static_cast<double>(k)));
}

} // namespace xsf
