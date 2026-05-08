#pragma once

#include "xsf/config.h"
#include "xsf/error.h"

/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds the functions
 * which make the recursion matrices.  This particular impl
 * returns three vectors, corresponding to the on and off diagonal
 * elements of the recursion matrix.
 *
 * Stuart Brorson, Summer 2025.
 *
 */

namespace xsf {
namespace mathieu {

    enum class Parity { Even, Odd };

    /* Get index of characteristic value in sorted array of eigenvalues
     * based on the order and parity. */

    template <Parity P, typename T>
    XSF_HOST_DEVICE T cv_index(T m) {
        static_assert(std::is_integral_v<T>, "m must be of integer type");

        if constexpr (P == Parity::Even) {
            return m / 2;
        }
        return (m - 1) / 2;
    }

    template <Parity FuncParity, Parity OrderParity>
    XSF_HOST_DEVICE double d0(double q) {
        constexpr auto Even = Parity::Even;
        constexpr auto Odd = Parity::Odd;

        if constexpr (FuncParity == Even && OrderParity == Even) {
            return 0.0;
        }
        if constexpr (FuncParity == Even && OrderParity == Odd) {
            return 1.0 + q;
        }
        if constexpr (FuncParity == Odd && OrderParity == Even) {
            return 4.0;
        }
        // FuncParity == Odd && OrderParity == Odd
        return 1.0 - q;
    }

    template <Parity FuncParity, Parity OrderParity>
    XSF_HOST_DEVICE double e0(double q) {
        constexpr auto Even = Parity::Even;

        if constexpr (FuncParity == Even && OrderParity == Even) {
            return M_SQRT2 * q;
        }
        return q;
    }

    template <Parity FuncParity, Parity OrderParity, typename U>
    XSF_HOST_DEVICE double sqrt_di(U j) {
        static_assert(std::is_integral_v<U>, "j must be of integer type");
        constexpr auto Even = Parity::Even;
        constexpr auto Odd = Parity::Odd;
        if constexpr (OrderParity == Odd) {
            return 2.0 * static_cast<double>(j) + 1.0;
        }
        if constexpr (FuncParity == Even) {
            return 2.0 * static_cast<double>(j);
        }
        // FuncParity == Odd && OrderParity == Even
        return 2.0 * (static_cast<double>(j) + 1.0);
    }

    /*-----------------------------------------------
      This creates the recurrence relation matrix for
      the even-even (ce_2n), even-odd (ce_2n+1),
      (se_2n) -- sometimes called se_2n+2, and odd-odd (se_2n+1)
      Mathieu functions, depending on the values of the template
      arguments FuncParity and OrderParity.

      Template arguments:

      FuncParity = parity of mathieu function, Parity::Even for the
      even Mathieu function and Parity::Odd for the odd
      Mathieu function.

      OrderParity = parity of the order m of a Mathieu function.

      Inputs:
      q = shape parameter.

      Outputs:
      D = Diagonal elements in recurrence matrix. mdspan view of a 1d array of
      size n.
      E = Off-diagonal elements in recurrence matrix. mdspan view of 1d array of
      size n - 1.
      -------------------------------------------------*/

    template <Parity FuncParity, Parity OrderParity, typename OutputMat>
    XSF_HOST_DEVICE void make_matrix(typename OutputMat::value_type q, OutputMat D, OutputMat E) {
        auto n = D.extent(0);
        if (n == 0) {
            return;
        }

        D(0) = d0<FuncParity, OrderParity>(q);

        for (decltype(n) j = 1; j < n; j++) {
            auto tj = sqrt_di<FuncParity, OrderParity>(j);
            D(j) = tj * tj;
        }

        if (n == 1) {
            return;
        }

        // Make off-diagonal entries
        E(0) = e0<FuncParity, OrderParity>(q);
        for (decltype(n) j = 1; j < (n - 1); j++) {
            E(j) = q;
        }
    }

    // Determine size of recursion matrix to use.
    XSF_HOST_DEVICE inline int get_partial_sum_N(int m, double q) {
        int N;

        // This is sort of ad-hoc ...
        if (q > 1.0) {
            double qq = std::log10(q); // I need to use size of q to compute N.
            N = m + 25 + 10 * qq;
        } else {
            N = m + 25;
        }
        return N;
    }

    template <Parity FuncParity, Parity OrderParity, typename InputMat>
    XSF_HOST_DEVICE void sum_fourier_series(
        InputMat X, double v, typename InputMat::value_type &out, typename InputMat::value_type &out_diff
    ) {
        auto constexpr Even = Parity::Even;
        auto constexpr Odd = Parity::Odd;
        auto N = X.extent(0);
        // Local scope variables used in summing the Fourier series.
        double tt, td, xep{0.0}, xem{0.0}, xedp{0.0}, xedm{0.0};

        // Sum from smallest to largest coeff.
        for (decltype(N) kp1 = N; kp1 > 0; kp1--) {
            auto k = kp1 - 1;
            auto r = sqrt_di<FuncParity, OrderParity>(k);
            auto phi = r * v;
            if constexpr (FuncParity == Even) {
                tt = X(k) * std::cos(phi);
                td = -r * X(k) * std::sin(phi);
            } else {
                tt = X(k) * std::sin(phi);
                td = r * X(k) * std::cos(phi);
            }
            if (tt < 0) {
                xem = xem + tt; // Neg running sum
            } else {
                xep = xep + tt; // Pos running sum
            }

            if (td < 0) {
                xedm = xedm + td;
            } else {
                xedp = xedp + td;
            }
        } // for (k=(N-1) ...

        // This makes sure the fcn has the right overall sign for q<0.
        // Someday combine this with the above sums into the same for loop.
        double x = 0.0;
        for (decltype(N) l = 0; l < N; l++) {
            x += X(l);
        }
        out = std::copysign(xep + xem, x);
        out_diff = std::copysign(xedp + xedm, x);
    }

} // namespace mathieu
} // namespace xsf
