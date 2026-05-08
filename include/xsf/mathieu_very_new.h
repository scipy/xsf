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

        if constexpr(P == Parity::Even) {
            return m / 2;
        }
        return (m - 1) / 2;
    }

    template <Parity FuncParity, Parity OrderParity, typename T>
    XSF_HOST_DEVICE T d0(T q) {
        constexpr auto Even = Parity::Even;
        constexpr auto Odd = Parity::Odd;

        if constexpr(FuncParity == Even && OrderParity == Even) {
            return T(0);
        }
        if constexpr(FuncParity == Even && OrderParity == Odd) {
            return T(1) + q;
        }
        if constexpr(FuncParity == Odd && OrderParity == Even) {
            return T(4);
        }
        // FuncParity == Odd && OrderParity == Odd
        return T(1) - q;
    }

    template <Parity FuncParity, Parity OrderParity, typename T>
    XSF_HOST_DEVICE T e0(T q) {
        constexpr auto Even = Parity::Even;

        if constexpr(FuncParity == Even && OrderParity == Even) {
            return T(M_SQRT2) * q;
        }
        return q;
    }

    template <Parity FuncParity, Parity OrderParity, typename T, typename U>
    XSF_HOST_DEVICE T sqrt_di(U j) {
        static_assert(std::is_integral_v<U>, "j must be of integer type");
        constexpr auto Even = Parity::Even;
        constexpr auto Odd = Parity::Odd;
        if constexpr(OrderParity == Odd) {
            return T(2) * static_cast<T>(j) + T(1);
        }
        if constexpr(FuncParity == Even) {
            return T(2) * static_cast<T>(j);
        }
        // FuncParity == Odd && OrderParity == Even
        return T(2) * (static_cast<T>(j) + T(1));
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

    template <Parity FuncParity, Parity OrderParity, typename T, typename OutputMat>
    XSF_HOST_DEVICE void make_matrix(T q, OutputMat D, OutputMat E) {
        auto n = D.extent(0);
        if (n == 0) {
            return;
        }

        D(0) = d0<FuncParity, OrderParity>(q);


        for (decltype(n) j = 1; j < n; j++) {
            auto tj = sqrt_di<FuncParity, OrderParity, T>(j);
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

} // namespace mathieu
} // namespace xsf
