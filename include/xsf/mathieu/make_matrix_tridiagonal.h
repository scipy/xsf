#ifndef MAKE_MATRIX_TRIDIAGONAL_H
#define MAKE_MATRIX_TRIDIAGONAL_H

#include "../config.h"
#include "../error.h"
#include "matrix_utils.h"

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

// #define SQRT2 1.414213562373095

namespace xsf {
namespace mathieu {

    /*-----------------------------------------------
      This creates the recurrence relation matrix for
      the even-even Mathieu fcns (ce_2n).

      Inputs:
      N = matrix size (related to max order desired).
      q = shape parameter.

      Outputs
      D = Diagonal elements in recurrence matrix (must be calloc'ed in caller).
      E = Off-diagonal elements in recurrence matrix (must be calloc'ed in caller).

      Return:
      return code = SF_ERROR_OK if OK.
      -------------------------------------------------*/
    int make_matrix_ee(int N, double q, double *D, double *E) {
        int j;

        // Make diagonal entries
        D[0] = 0.0;
        for (j = 1; j < N; j++) {
            D[j] = (2.0 * j) * (2.0 * j);
        }

        // Make off-diagonal entries
        E[0] = M_SQRT2 * q;
        for (j = 1; j < (N - 1); j++) {
            E[j] = q;
        }

        return SF_ERROR_OK;
    }

    /*-----------------------------------------------
      This creates the recurrence relation matrix for the
      even-odd Mathieu fcns (ce_2n+1).

      Inputs:
      N = matrix size (related to max order desired).
      q = shape parameter.

      Outputs:
      D = Diagonal elements in recurrence matrix (must be calloc'ed in caller).
      E = Off-diagonal elements in recurrence matrix (must be calloc'ed in caller).

      Return:
      return code = SF_ERROR_OK if OK.
      -------------------------------------------------*/
    int make_matrix_eo(int N, double q, double *D, double *E) {
        int j;

        // Make diagonal entries
        D[0] = 1.0 + q;
        for (j = 1; j < N; j++) {
            D[j] = (2.0 * j + 1.0) * (2.0 * j + 1.0);
        }

        // Make off-diagonal entries
        for (j = 0; j < (N - 1); j++) {
            E[j] = q;
        }

        return SF_ERROR_OK;
    }

    /*-----------------------------------------------
      This creates the recurrence relation matrix for
      the odd-even Mathieu fcns (se_2n) -- sometimes called
      se_2n+2.

      Inputs:
      N = matrix size (related to max order desired).
      q = shape parameter.

      Outputs:
      D = Diagonal elements in recurrence matrix (must be calloc'ed in caller).
      E = Off-diagonal elements in recurrence matrix (must be calloc'ed in caller).

      Return:
      return code = SF_ERROR_OK if OK.
      -------------------------------------------------*/
    int make_matrix_oe(int N, double q, double *D, double *E) {
        int j;

        // Make diagonal entries
        for (j = 0; j < N; j++) {
            D[j] = (2.0 * (j + 1)) * (2.0 * (j + 1));
        }

        // Make off-diagonal entries
        for (j = 0; j < (N - 1); j++) {
            E[j] = q;
        }

        return SF_ERROR_OK;
    }

    /*-----------------------------------------------
      This creates the recurrence relation matrix for
      the odd-odd Mathieu fcns (se_2n+1).

      Inputs:
      N = matrix size (related to max order desired).
      q = shape parameter.

      Outputs:
      D = Diagonal elements in recurrence matrix (must be calloc'ed in caller).
      E = Off-diagonal elements in recurrence matrix (must be calloc'ed in caller).

      Return:
      return code = SF_ERROR_OK if OK.
      -------------------------------------------------*/
    int make_matrix_oo(int N, double q, double *D, double *E) {
        int j;

        // Make diagonal entries
        D[0] = 1 - q;
        for (j = 1; j < N; j++) {
            D[j] = (2.0 * j + 1.0) * (2.0 * j + 1.0);
        }

        // Make off-diagonal entries
        for (j = 0; j < (N - 1); j++) {
            E[j] = q;
        }

        return SF_ERROR_OK;
    }

} // namespace mathieu
} // namespace xsf

#endif // #ifndef MAKE_MATRIX_TRIDIAGONAL_H
