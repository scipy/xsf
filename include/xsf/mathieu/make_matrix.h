#ifndef MAKE_MATRIX_H
#define MAKE_MATRIX_H

#include "../config.h"
#include "../error.h"
#include "matrix_utils.h"

/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds the functions
 * which make the recursion matrices.
 *
 * Stuart Brorson, Summer 2025.
 *
 */

#define SQRT2 1.414213562373095

namespace xsf {
namespace mathieu {

    /*-----------------------------------------------
    This creates the recurrence relation matrix for
    the even-even Mathieu fcns (ce_2n).
    Inputs:
    N = matrix size (related to max order desired).
    q = shape parameter.
    Output:
    A = recurrence matrix (must be calloc'ed in caller).
    Return:
    return code = 0 if OK.
    -------------------------------------------------*/
    int make_matrix_ee(int N, double q, double *A) {
        int j;
        int i;

        // Symmetrize matrix here, then fix in caller.
        i = MATRIX_IDX(N, 0, 1);
        A[i] = SQRT2 * q;
        i = MATRIX_IDX(N, 1, 0);
        A[i] = SQRT2 * q;
        i = MATRIX_IDX(N, 1, 1);
        A[i] = 4.0;
        i = MATRIX_IDX(N, 1, 2);
        A[i] = q;

        for (j = 2; j <= N - 2; j++) {
            i = MATRIX_IDX(N, j, j - 1);
            A[i] = q;
            i = MATRIX_IDX(N, j, j);
            A[i] = (2.0 * j) * (2.0 * j);
            i = MATRIX_IDX(N, j, j + 1);
            A[i] = q;
        }

        i = MATRIX_IDX(N, N - 1, N - 2);
        A[i] = q;
        i = MATRIX_IDX(N, N - 1, N - 1);
        A[i] = (2.0 * (N - 1)) * (2.0 * (N - 1));

        return 0;
    }

    /*-----------------------------------------------
    This creates the recurrence relation matrix for
    the even-odd Mathieu fcns (ce_2n+1).
    Inputs:
    N = matrix size (related to max order desired).
    q = shape parameter.
    Output:
    A = recurrence matrix (calloc in caller).
    Return:
    return code = 0 if OK.
    -------------------------------------------------*/
    int make_matrix_eo(int N, double q, double *A) {
        int j;
        int i;

        i = MATRIX_IDX(N, 0, 0);
        A[i] = 1.0 + q;
        i = MATRIX_IDX(N, 0, 1);
        A[i] = q;
        i = MATRIX_IDX(N, 1, 0);
        A[i] = q;
        i = MATRIX_IDX(N, 1, 1);
        A[i] = 9.0;
        i = MATRIX_IDX(N, 1, 2);
        A[i] = q;

        for (j = 2; j <= N - 2; j++) {
            i = MATRIX_IDX(N, j, j - 1);
            A[i] = q;
            i = MATRIX_IDX(N, j, j);
            A[i] = (2.0 * j + 1.0) * (2.0 * j + 1.0);
            i = MATRIX_IDX(N, j, j + 1);
            A[i] = q;
        }

        i = MATRIX_IDX(N, N - 1, N - 2);
        A[i] = q;
        i = MATRIX_IDX(N, N - 1, N - 1);
        A[i] = (2.0 * (N - 1) + 1.0) * (2.0 * (N - 1) + 1.0);

        return 0;
    }

    /*-----------------------------------------------
    This creates the recurrence relation matrix for
    the odd-even Mathieu fcns (se_2n) -- sometimes called
    se_2n+2.
    Inputs:
    N = matrix size (related to max order desired).
    q = shape parameter.
    Output:
    A = recurrence matrix (calloc in caller).
    Return:
    return code = 0 if OK.
    -------------------------------------------------*/
    int make_matrix_oe(int N, double q, double *A) {
        int j;
        int i;

        i = MATRIX_IDX(N, 0, 0);
        A[i] = 4.0;
        i = MATRIX_IDX(N, 0, 1);
        A[i] = q;
        i = MATRIX_IDX(N, 1, 0);
        A[i] = q;
        i = MATRIX_IDX(N, 1, 1);
        A[i] = 16.0;
        i = MATRIX_IDX(N, 1, 2);
        A[i] = q;

        for (j = 2; j <= N - 2; j++) {
            i = MATRIX_IDX(N, j, j - 1);
            A[i] = q;
            i = MATRIX_IDX(N, j, j);
            A[i] = (2.0 * (j + 1)) * (2.0 * (j + 1));
            i = MATRIX_IDX(N, j, j + 1);
            A[i] = q;
        }

        i = MATRIX_IDX(N, N - 1, N - 2);
        A[i] = q;
        i = MATRIX_IDX(N, N - 1, N - 1);
        A[i] = (2.0 * N) * (2.0 * N);

        return 0;
    }

    /*-----------------------------------------------
    This creates the recurrence relation matrix for
    the odd-odd Mathieu fcns (se_2n+1).
    Inputs:
    N = matrix size (related to max order desired).
    q = shape parameter.
    Output:
    A = recurrence matrix (calloc in caller).
    Return:
    return code = 0 if OK.
    -------------------------------------------------*/
    int make_matrix_oo(int N, double q, double *A) {
        int j;
        int i;

        i = MATRIX_IDX(N, 0, 0);
        A[i] = 1.0 - q;
        i = MATRIX_IDX(N, 0, 1);
        A[i] = q;
        i = MATRIX_IDX(N, 1, 0);
        A[i] = q;
        i = MATRIX_IDX(N, 1, 1);
        A[i] = 9.0;
        i = MATRIX_IDX(N, 1, 2);
        A[i] = q;

        for (j = 2; j <= N - 2; j++) {
            i = MATRIX_IDX(N, j, j - 1);
            A[i] = q;
            i = MATRIX_IDX(N, j, j);
            A[i] = (2.0 * j + 1.0) * (2.0 * j + 1.0);
            i = MATRIX_IDX(N, j, j + 1);
            A[i] = q;
        }

        i = MATRIX_IDX(N, N - 1, N - 2);
        A[i] = q;
        i = MATRIX_IDX(N, N - 1, N - 1);
        A[i] = (2.0 * N - 1.0) * (2.0 * N - 1.0);

        return 0;
    }

} // namespace mathieu
} // namespace xsf

#endif // #ifndef MAKE_MATRIX_H
