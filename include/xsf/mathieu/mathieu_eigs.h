#ifndef MATHIEU_EIGS_H
#define MATHIEU_EIGS_H

#include "../config.h"
#include "../error.h"
#include "make_matrix.h"
#include "matrix_utils.h"
#include <vector>

/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds the functions
 * which return the Mathieu eigenvalues (characteristic values) a and
 * b as a function of parameter q.
 *
 * Stuart Brorson, summer 2025.
 *
 */

/* DSYEVD_ prototype */
#ifdef __cplusplus
extern "C" {
#endif
void dsyevd_(
    char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *iwork, int *liwork,
    int *info
);
#ifdef __cplusplus
}
#endif

namespace xsf {
namespace mathieu {

    //------------------------------------------------------
    // This is the Mathieu characteristic value (eigenvalue)
    // a for even fcns.
    int mathieu_a(int m, double q, double *a) {
        // printf("--> mathieu_a, m = %d, q = %e\n", m, q);

        int N = m + 25; // Sets size of recursion matrix
        int retcode = SF_ERROR_OK;

        if (m > 500) {
            // Don't support absurdly larger orders for now.
            *a = std::numeric_limits<double>::quiet_NaN();
            return SF_ERROR_DOMAIN;
        }

        // Allocate recursion matrix
        std::vector<double> A(N * N);

        // Allocate vector for eigenvalues
        std::vector<double> ww(N);

        // Do EVD
        if (m % 2 == 0) {
            // Even order m
            retcode = make_matrix_ee(N, q, A.data());
            if (retcode != SF_ERROR_OK) {
                *a = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_OTHER; // Not sure what went wrong.
            }

            char V = 'N';
            char U = 'U';

            /* Query and allocate the optimal workspace */
            int lwork = -1;
            int liwork = -1;
            double work_query;
            int iwork_query;
            dsyevd_(&V, &U, &N, A.data(), &N, ww.data(), &work_query, &lwork, &iwork_query, &liwork, &retcode);
            lwork = (int)work_query;
            liwork = iwork_query;

            /* Allocate worksapce */
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);

            /* Solve eigenproblem */
            dsyevd_(&V, &U, &N, A.data(), &N, ww.data(), work.data(), &lwork, iwork.data(), &liwork, &retcode);

            // Check if dsyevd was successful
            if (retcode != 0) {
                *a = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_NO_RESULT;
            }

            // Sort ww vector from lowest to highest
            // quickSort(ww, 0, N-1);
            // print_matrix(ww, N, 1);

            // Now figure out which one to return.
            int idx = m / 2;
            *a = ww[idx];

        } else {
            // Odd order m
            retcode = make_matrix_eo(N, q, A.data());
            if (retcode != SF_ERROR_OK) {
                *a = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_OTHER;
            }

            char V = 'V';
            char U = 'U';

            /* Query and allocate the optimal workspace */
            int lwork = -1;
            int liwork = -1;
            double work_query;
            int iwork_query;
            dsyevd_(&V, &U, &N, A.data(), &N, ww.data(), &work_query, &lwork, &iwork_query, &liwork, &retcode);
            lwork = (int)work_query;
            liwork = iwork_query;

            /* Allocate worksapce */
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);

            /* Solve eigenproblem */
            dsyevd_(&V, &U, &N, A.data(), &N, ww.data(), work.data(), &lwork, iwork.data(), &liwork, &retcode);

            // Check if dsyevd was successful
            if (retcode != 0) {
                *a = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_NO_RESULT;
            }

            // Sort ww vector from lowest to highest
            // quickSort(ww, 0, N-1);
            // print_matrix(ww, N, 1);

            // Now figure out which one to return.
            int idx = (m - 1) / 2;
            *a = ww[idx];
        }

        // printf("<-- mathieu_a\n");
        return retcode;
    }

    //------------------------------------------------------
    int mathieu_b(int m, double q, double *b) {
        // This computes the Mathieu characteristic value (eigenvalue)
        // for odd fcns.
        // printf("--> mathieu_b, m = %d, q = %e\n", m, q);
        int N = m + 25; // Sets size of recursion matrix
        int retcode = SF_ERROR_OK;

        if (m > 500) {
            // Don't support absurdly larger orders for now.
            *b = std::numeric_limits<double>::quiet_NaN();
            return SF_ERROR_DOMAIN;
        }

        // Allocate recursion matrix
        std::vector<double> B(N * N);

        // Allocate vector for eigenvalues
        std::vector<double> ww(N);

        // Do EVD
        if (m % 2 == 0) {
            // Even order m
            retcode = make_matrix_oe(N, q, B.data());
            if (retcode != SF_ERROR_OK) {
                *b = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_OTHER;
            }

            char V = 'V';
            char U = 'U';

            /* Query and allocate the optimal workspace */
            int lwork = -1;
            int liwork = -1;
            double work_query;
            int iwork_query;
            dsyevd_(&V, &U, &N, B.data(), &N, ww.data(), &work_query, &lwork, &iwork_query, &liwork, &retcode);
            lwork = (int)work_query;
            liwork = iwork_query;

            /* Allocate worksapce */
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);

            /* Solve eigenproblem */
            dsyevd_(&V, &U, &N, B.data(), &N, ww.data(), work.data(), &lwork, iwork.data(), &liwork, &retcode);

            if (retcode != 0) {
                *b = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_NO_RESULT;
            }

            // Sort ww vector from lowest to highest
            // quickSort(ww, 0, N-1);
            // print_matrix(ww, N, 1);

            // Now figure out which one to return.
            int idx = (m - 2) / 2;
            *b = ww[idx];

        } else {
            // Odd order m
            retcode = make_matrix_oo(N, q, B.data());
            if (retcode != SF_ERROR_OK) {
                *b = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_OTHER;
            }

            char V = 'V';
            char U = 'U';

            /* Query and allocate the optimal workspace */
            int lwork = -1;
            int liwork = -1;
            double work_query;
            int iwork_query;
            dsyevd_(&V, &U, &N, B.data(), &N, ww.data(), &work_query, &lwork, &iwork_query, &liwork, &retcode);
            lwork = (int)work_query;
            liwork = iwork_query;

            /* Allocate worksapce */
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);

            /* Solve eigenproblem */
            dsyevd_(&V, &U, &N, B.data(), &N, ww.data(), work.data(), &lwork, iwork.data(), &liwork, &retcode);

            if (retcode != 0) {
                *b = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_NO_RESULT;
            }

            // Sort ww vector from lowest to highest
            // quickSort(ww, 0, N-1);
            // print_matrix(ww, N, 1);

            // Now figure out which one to return.
            int idx = (m - 1) / 2;
            *b = ww[idx];
        }

        // printf("<-- mathieu_b\n");
        return retcode;
    }

} // namespace mathieu
} // namespace xsf

#endif // #ifndef MATHIEU_EIGS_H
