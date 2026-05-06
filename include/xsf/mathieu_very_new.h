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

    /*-----------------------------------------------
      This creates the recurrence relation matrix for
      the even-even Mathieu fcns (ce_2n).

      Inputs:
      q = shape parameter.

      Outputs:
      D = Diagonal elements in recurrence matrix. mdspan view of a 1d array of
      size n.
      E = Off-diagonal elements in recurrence matrix. mdspan view of 1d array of
      size n - 1.
      -------------------------------------------------*/

    template <typename T, typename OutputMat>
    XSF_HOST_DEVICE inline void make_matrix_ee(T q, OutputMat D, OutputMat E) {
        auto n = D.extent(0);

        // Make diagonal entries
        D(0) = T(0);
        for (decltype(n) j = 1; j < n; j++) {
            auto tj = T(2) * static_cast<T>(j);
            D(j) = tj * tj;
        }

        // Make off-diagonal entries
        E(0) = T(M_SQRT2) * q;
        for (decltype(n) j = 1; j < (n - 1); j++) {
            E(j) = q;
        }
    }

    /*-----------------------------------------------
      This creates the recurrence relation matrix for the
      even-odd Mathieu fcns (ce_2n+1).

      Inputs:
      q = shape parameter.

      Outputs:
      D = Diagonal elements in recurrence matrix. mdspan view of a 1d array of
      size n.
      E = Off-diagonal elements in recurrence matrix. mdspan view of 1d array of
      size n - 1.
      -------------------------------------------------*/
    template <typename T, typename OutputMat>
    XSF_HOST_DEVICE inline void make_matrix_eo(T q, OutputMat D, OutputMat E) {
        auto n = D.extent(0);

        // Make diagonal entries
        D(0) = T(1) + q;
        for (decltype(n) j = 1; j < n; j++) {
            auto tj = T(2) * static_cast<T>(j) + 1.0;
            D(j) = tj * tj;
        }

        // Make off-diagonal entries
        for (decltype(n) j = 0; j < (n - 1); j++) {
            E(j) = q;
        }
    }

    /*-----------------------------------------------
      This creates the recurrence relation matrix for
      the odd-even Mathieu fcns (se_2n) -- sometimes called
      se_2n+2.

      Inputs:
      q = shape parameter.

      Outputs:
      D = Diagonal elements in recurrence matrix. mdspan view of a 1d array of size
      n.
      E = Off-diagonal elements in recurrence matrix. mdspan view of 1d array of size
      n - 1.
      -------------------------------------------------*/
    template <typename T, typename OutputMat>
    XSF_HOST_DEVICE inline void make_matrix_oe(T q, OutputMat D, OutputMat E) {
        // Coeffs for se_2n+2
        auto n = D.extent(0);

        // Make diagonal entries
        for (decltype(n) j = 0; j < n; j++) {
            auto tj = T(2) * (static_cast<T>(j) + T(1));
            D(j) = tj * tj;
        }

        // Make off-diagonal entries
        for (decltype(n) j = 0; j < (n - 1); j++) {
            E(j) = q;
        }
    }

    /*-----------------------------------------------
      This creates the recurrence relation matrix for
      the odd-odd Mathieu fcns (se_2n+1).

      Inputs:
      q = shape parameter.

      D = Diagonal elements in recurrence matrix. mdspan view of a 1d array of size
      n.
      E = Off-diagonal elements in recurrence matrix. mdspan view of 1d array of size
      n - 1.
      -------------------------------------------------*/
    template <typename T, typename OutputMat>
    XSF_HOST_DEVICE inline void make_matrix_oo(T q, OutputMat D, OutputMat E) {
        // Coeffs for se_2n+1
        auto n = D.extent(0);

        // Make diagonal entries
        D(0) = T(1) - q;
        for (decltype(n) j = 1; j < n; j++) {
            auto tj = T(2) * static_cast<T>(j) + T(1);
            D(j) = tj * tj;
        }

        // Make off-diagonal entries
        for (decltype(n) j = 0; j < (n - 1); j++) {
            E(j) = q;
        }
    }

    int get_partial_sum_N(int m, double q) {
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
