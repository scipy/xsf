#include "../testing_utils.h"

#include <xsf/mathieu.h>

#include <cstdio>
#include <stdlib.h>

/*
 *
 * The goal of these tests are to verify that my C/C++
 * impl of the Mathieu fcns has been carried over from
 * my Matlab impl correctly.  Therefore, main just calls
 * a bunch of golden value tests.  The GVs were generated
 * using the Matlab impl.  I did fairly extensive
 * validation of the Matlab impl.  Therefore, if the tests
 * here show the C impl matches the Matlab impl, then that
 * should serve as verification of the C impl's correctness.
 *
 * A secondary goal is to show how to call the various fcns
 * in my API.
 *
 */

namespace xsf {
namespace mathieu {

    //-------------------------------------------------------------
    extern "C" int main() {
        int N = 6;
        int pass = 0;
        int fail = 0;

        //*******************************************************
        // First print out the recursion matrices.  These are
        // private -- not public fcns.  But I want to see they
        // are correct.
        //*******************************************************
        {
            double *A = (double *)calloc(N * N, sizeof(double));
            double q = 2.0;

            make_matrix_ee(N, q, A);
            print_matrix(A, N, N);
            printf("----------------------------------------------\n");

            make_matrix_eo(N, q, A);
            print_matrix(A, N, N);
            printf("----------------------------------------------\n");

            make_matrix_oe(N, q, A);
            print_matrix(A, N, N);
            printf("----------------------------------------------\n");

            make_matrix_oo(N, q, A);
            print_matrix(A, N, N);
            printf("----------------------------------------------\n");

            free(A);
        }

        //*******************************************************
        // Try computing the eigenvalues.  These are public fcns.
        //*******************************************************
        ///*
        printf("==============================================\n");
        printf("Test a eigenvalues\n");
        double a;
        {
            double q = 0.001;
            double tol = 1e-13;

            // Golden values from Matlab.
            double a_true[6] = {-4.999999453125127e-07, 1.000999874984374,  4.000000416666611,
                                9.000000062515628,      16.000000033333333, 25.000000020833337};
            printf("q = %f\n", q);
            for (int m = 0; m < N; m++) {
                mathieu_a(m, q, &a);
                printf("Order m = %d, eigenvalue a = %8.5f\n", m, a);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(a - a_true[m]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            double a_true[6] = {-0.455138604107414, 1.859108072514364,  4.371300982735085,
                                9.078368847203102,  16.033832340359510, 25.020854345448583};
            printf("q = %f\n", q);
            for (int m = 0; m < N; m++) {
                mathieu_a(m, q, &a);
                printf("Order m = %d, eigenvalue a = %8.5f\n", m, a);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(a - a_true[m]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 100.0;
            double tol = 1e-13;
            double a_true[6] = {-1.802532491522514e+02, -1.412800568086200e+02, -1.033705070649672e+02,
                                -66.574389967485317,    -30.950104000623185,    3.432562047992776};
            q = 100.0;
            printf("q = %f\n", q);
            for (int m = 0; m < N; m++) {
                mathieu_a(m, q, &a);
                printf("Order m = %d, eigenvalue a = %8.5f\n", m, a);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(a - a_true[m]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }
        //*/

        ///*
        printf("==============================================\n");
        double b;
        printf("Test b eigenvalues\n");
        {
            double q = 0.001;
            double tol = 1e-13;
            double b_true[6] = {0.998999875015624,  3.999999916666667,  9.000000062484375,
                                16.000000033333333, 25.000000020833333, 36.000000014285725};
            printf("q = %f\n", q);
            for (int m = 1; m < N; m++) {
                mathieu_b(m, q, &b);
                printf("Order m = %d, eigenvalue b = %8.5f\n", m, b);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(b - b_true[m - 1]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            double b_true[6] = {-0.110248816992095, 3.917024772998471,  9.047739259809376,
                                16.032970081405796, 25.020840823289770, 36.014289910628236};
            printf("q = %f\n", q);
            for (int m = 1; m < N; m++) {
                mathieu_b(m, q, &b);
                printf("Order m = %d, eigenvalue b = %8.5f\n", m, b);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(b - b_true[m - 1]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 100.0;
            double tol = 1e-13;
            double b_true[6] = {-1.802532491522514e+02, -1.412800568086194e+02, -1.033705070649312e+02,
                                -66.574389965842357,    -30.950103947238059,    3.432563359324842};
            printf("q = %f\n", q);
            for (int m = 1; m < N; m++) {
                mathieu_b(m, q, &b);
                printf("Order m = %d, eigenvalue b = %8.5f\n", m, b);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(b - b_true[m - 1]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }
        //*/

        //*******************************************************
        // Compute the Fourier coeffs -- these are public fcns.
        //*******************************************************
        printf("==============================================\n");
        printf("Test A Fourier coeffs\n");
        int retcode;
        int m;
        double *AA = (double *)calloc(N, sizeof(double));
        N = 5;
        {
            double q = 1.0;
            double tol = 1e-13;
            double AA_true[5] = {
                6.72989672316633203e-01, -3.06303580036475842e-01, 1.86455593633495648e-02, -5.11683638542392003e-04,
                7.93860116701059744e-06
            };

            m = 0;
            printf("ee, N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_ee(N, q, m, AA);
            if (retcode == 0) {
                pass++;
                printf("Function call succeeded!\n");
            } else {
                printf("Call to mathieu_coeffs_ee returned failure!\n");
                fail++;
            }

            printf("Checking values ... \n");
            for (int j = 0; j < N; j++) {
                printf("AA[%d] = % 9.6e ...", j, AA[j]);
                // I test the relative err here since the coeffs vary widely
                if (abs((AA[j] - AA_true[j]) / AA[j]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! AA_true = % 9.6e\n", AA_true[j]);
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            m = 1;
            double AA_true[5] = {
                9.90202059407947477e-01, -1.39511476750210750e-01, 6.03431870922933183e-03, -1.28040356069675752e-04,
                1.61787860802725383e-06
            };
            printf("eo, N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_eo(N, q, m, AA);
            if (retcode == 0) {
                pass++;
                printf("Function call succeeded!\n");
            } else {
                printf("Call to mathieu_coeffs_eo returned failure!\n");
                fail++;
            }

            printf("Checking values ... \n");
            for (int j = 0; j < N; j++) {
                printf("AA[%d] = % 9.6e ...", j, AA[j]);
                if (abs((AA[j] - AA_true[j]) / AA[j]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! AA_true = % 9.6e\n", AA_true[j]);
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            m = 2;
            double AA_true[5] = {
                2.16927946751974854e-01, 9.48257346823881631e-01, -8.17670087238141774e-02, 2.58658716581847224e-03,
                -4.33782257276887011e-05
            };
            printf("N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_ee(N, q, m, AA);
            if (retcode == 0) {
                pass++;
                printf("Function call succeeded!\n");
            } else {
                printf("Call to mathieu_coeffs_ee returned failure!\n");
                fail++;
            }

            printf("Checking values ... \n");
            for (int j = 0; j < N; j++) {
                printf("AA[%d] = % 9.6e ...", j, AA[j]);
                if (abs((AA[j] - AA_true[j]) / AA[j]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! AA_true = % 9.6e\n", AA_true[j]);
                }
            }
        }

        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            m = 3;
            double AA_true[5] = {
                1.39615654604830663e-01, 9.88251100137287120e-01, -6.21675551357310924e-02, 1.55778240472689887e-03,
                -2.16594420865885830e-05
            };

            printf("N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_eo(N, q, m, AA);
            if (retcode == 0) {
                pass++;
                printf("Function call succeeded!\n");
            } else {
                printf("Call to mathieu_coeffs_eo returned failure!\n");
                fail++;
            }

            printf("Checking values ... \n");
            for (int j = 0; j < N; j++) {
                printf("AA[%d] = % 9.6e ...", j, AA[j]);
                if (abs((AA[j] - AA_true[j]) / AA[j]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! AA_true = % 9.6e\n", AA_true[j]);
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            printf("This should fail -- calling coeffs_eo with even m\n");
            double q = 1.0;
            m = 4;
            printf("N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_eo(N, q, m, AA);
            if (retcode != 0) {
                pass++;
                printf("Function call correctly returned error!\n");
            } else {
                printf("Bad call to mathieu_coeffs_eo returned success -- bad!\n");
                fail++;
            }
        }

        //=====================================================
        printf("==============================================\n");
        printf("Test B Fourier coeffs\n");
        {
            printf("This should fail -- calling coeffs_oe with m=0\n");
            double q = 1.0;
            m = 0;
            printf("N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_oe(N, q, m, AA);
            if (retcode != 0) {
                pass++;
                printf("Function call correctly returned error!\n");
            } else {
                printf("Bad call to mathieu_coeffs_eo returned success -- bad!\n");
                fail++;
            }
        }

        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            m = 1;
            double AA_true[5] = {
                9.93967961398935729e-01, -1.09583791872267272e-01, 4.36764886689439396e-03, -8.89579207045385789e-05,
                1.09675314774650835e-06
            };
            printf("N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_oo(N, q, m, AA);
            if (retcode == 0) {
                pass++;
                printf("Function call succeeded!\n");
            } else {
                printf("Call to mathieu_coeffs_oo returned failure!\n");
                fail++;
            }

            printf("Checking values ... \n");
            for (int j = 0; j < N; j++) {
                printf("AA[%d] = % 9.6e ...", j, AA[j]);
                if (abs((AA[j] - AA_true[j]) / AA[j]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! AA_true = % 9.6e\n", AA_true[j]);
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            m = 2;
            double AA_true[5] = {
                9.96571915618007398e-01,  -8.26907809217511391e-02, 2.57874176092240956e-03,
                -4.29271107568478840e-05, 4.46771248032548085e-07,
            };

            printf("N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_oe(N, q, m, AA);
            if (retcode == 0) {
                pass++;
                printf("Function call succeeded!\n");
            } else {
                printf("Call to mathieu_coeffs_oe returned failure!\n");
                fail++;
            }

            printf("Checking values ... \n");
            for (int j = 0; j < N; j++) {
                printf("AA[%d] = % 9.6e ...", j, AA[j]);
                if (abs((AA[j] - AA_true[j]) / AA[j]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! AA_true = % 9.6e\n", AA_true[j]);
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            m = 3;
            double AA_true[5] = {
                1.09642473301236554e-01, 9.92016510230659843e-01, -6.22843393799822065e-02, 1.55951158907772385e-03,
                -2.16742541934712708e-05
            };

            printf("N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_oo(N, q, m, AA);
            if (retcode == 0) {
                pass++;
                printf("Function call succeeded!\n");
            } else {
                printf("Call to mathieu_coeffs_oo returned failure!\n");
                fail++;
            }

            printf("Checking values ... \n");
            for (int j = 0; j < N; j++) {
                printf("AA[%d] = % 9.6e ...", j, AA[j]);
                if (abs((AA[j] - AA_true[j]) / AA[j]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! AA_true = % 9.6e\n", AA_true[j]);
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 1.0;
            double tol = 1e-13;
            m = 4;
            double AA_true[5] = {
                8.27162782989298156e-02, 9.95322502016357302e-01, -4.99004143812368725e-02, 1.04056488398845915e-03,
                -1.23925412747991493e-05
            };

            printf("N = %d, q = %f, m = %d\n", N, q, m);
            retcode = mathieu_coeffs_oe(N, q, m, AA);
            if (retcode == 0) {
                pass++;
                printf("Function call succeeded!\n");
            } else {
                printf("Call to mathieu_coeffs_oe returned failure!\n");
                fail++;
            }

            printf("Checking values ... \n");
            for (int j = 0; j < N; j++) {
                printf("AA[%d] = % 9.6e ...", j, AA[j]);
                if (abs((AA[j] - AA_true[j]) / AA[j]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! AA_true = % 9.6e\n", AA_true[j]);
                }
            }
        }

        //*******************************************************
        // Compute the Mathieu fcns & compare the results to
        // those from the Matlab impls.  These are
        // the publically-accessible Mathieu fcn impls.
        //*******************************************************
        printf("==============================================\n");
        printf("Test Mathieu ce\n");
        double ce;
        double ced;
        N = 6;
        {
            double v = 1.0;
            double q = 0.001;
            double tol = 1e-13;
            // Golden values obtained from Matlab for m = 0, 1, 2, 3, 4, 5.
            double ce_true[6] = {0.707253852647939,  0.540426067653929,  -0.415842336334926,
                                 -0.989942668408815, -0.653726300137563, 0.283568898267299};
            double ced_true[6] = {6.430371575928823e-04, -0.841418026637729, -1.818846996799653,
                                  -0.423764888085533,    3.026974584470965,  4.794786515927561};

            printf("q = %f\n", q);
            for (int m = 0; m < N; m++) {
                retcode = mathieu_ce(m, q, v, &ce, &ced);
                printf("Order m = %d, v = %f, ce = %8.5f, ced = %8.5f\n", m, v, ce, ced);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(ce - ce_true[m]) < tol && abs(ced - ced_true[m]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double v = 1.0;
            double q = 100;
            double tol = 1e-13;
            double ce_true[6] = {0.086725459442521, 0.315928563892497, 0.737195938153019,
                                 1.233155548880747, 1.476584096867507, 1.114354874812387};
            double ced_true[6] = {
                0.924256230828739, 2.780839256567336, 4.896513216077115, 4.882886873958746,
                0.264149721572143, -7.943564372635874

            };
            printf("q = %f\n", q);
            for (int m = 0; m < N; m++) {
                retcode = mathieu_ce(m, q, v, &ce, &ced);
                printf("Order m = %d, v = %f, ce = %8.5f, ced = %8.5f\n", m, v, ce, ced);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(ce - ce_true[m]) < tol && abs(ced - ced_true[m]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }

        printf("==============================================\n");
        printf("Test Mathieu se\n");
        double se;
        double sed;
        N = 6;
        {
            double v = 1.0;
            double q = 0.001;
            double tol = 1e-13;
            // Golden values obtained from Matlab for m = 1, 2, 3, 4, 5, 6.
            double se_true[6] = {0.841453335445950,  0.909360489814828,  0.141285111200371,
                                 -0.756712745143668, -0.958942823901047, -0.279488670902052};
            double sed_true[6] = {0.540673509812678,  -0.832075773996413, -2.969998567644164,
                                  -2.614931881211228, 1.417905406870398,  5.760932545758129};

            printf("q = %f\n", q);
            for (int m = 1; m <= N; m++) {
                retcode = mathieu_se(m, q, v, &se, &sed);
                printf("Order m = %d, v = %f, se = %8.5f, sed = %8.5f\n", m, v, se, sed);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(se - se_true[m - 1]) < tol && abs(sed - sed_true[m - 1]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }
        printf("----------------------------------------------\n");
        {
            double v = 1.0;
            double q = 100;
            double tol = 1e-13;
            // Golden values obtained from Matlab for m = 1, 2, 3, 4, 5, 6.
            double se_true[6] = {0.086725459442519, 0.315928563892424, 0.737195938151456,
                                 1.233155548910162, 1.476584099964093, 1.114354958780249};
            double sed_true[6] = {0.924256230828753, 2.780839256567937, 4.896513216103436,
                                  4.882886874949449, 0.264149739778980, -7.943564564030449};

            printf("q = %f\n", q);
            for (int m = 1; m <= N; m++) {
                retcode = mathieu_se(m, q, v, &se, &sed);
                printf("Order m = %d, v = %f, se = %8.5f, sed = %8.5f\n", m, v, se, sed);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(se - se_true[m - 1]) < tol && abs(sed - sed_true[m - 1]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }

        //*******************************************************
        // Check that my Bessel fcns return the right thing
        //*******************************************************
        printf("==============================================\n");
        printf("Test my BesselJ implementation ... \n");
        double J;
        double Y;
        {
            // Golden values obtained from Matlab for z = 1 and m = -3, -2, -1, 0, 1, 2, 3..
            double z = 1.0;
            double tol = 1e-13;
            double J_true[7] = {-0.019563353982668, 0.114903484931901, -0.440050585744934, 0.765197686557967,
                                0.440050585744934,  0.114903484931901, 0.019563353982668};

            for (int m = -3; m <= 3; m++) {
                J = besselj(m, z);
                printf("Order m = %d, z = %f, J = %20.17f, J_true = %20.17e\n", m, z, J, J_true[m + 3]);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(J - J_true[m + 3]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }

        printf("------------------------------------------------\n");
        printf("Test my BesselY implementation, z = 1 ... \n");
        {
            // Golden values obtained from Matlab for z = 1 and m = -3, -2, -1, 0, 1, 2, 3..
            double z = 1.0;
            double tol = 1e-13;
            double Y_true[7] = {5.821517605964731,  -1.650682606816255, 0.781212821300289, 0.088256964215677,
                                -0.781212821300289, -1.650682606816255, -5.821517605964731};

            for (int m = -3; m <= 3; m++) {
                Y = bessely(m, z);
                printf("Order m = %d, z = %f, Y = %20.17e, Y_true = %20.17e\n", m, z, Y, Y_true[m + 3]);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(Y - Y_true[m + 3]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }

        printf("------------------------------------------------\n");
        printf("Test my BesselY implementation, z = .01 ... \n");
        {
            // Golden values obtained from Matlab for z = 1 and m = -3, -2, -1, 0, 1, 2, 3..
            double z = 0.01;
            double tol = 1e-13;
            double Y_true[7] = {5.093021841713738e+06, -1.273271380077505e+04, 63.678596282060660,
                                -3.005455637083646,    -63.678596282060660,    -1.273271380077505e+04,
                                -5.093021841713738e+06};

            for (int m = -3; m <= 3; m++) {
                Y = bessely(m, z);
                printf("Order m = %d, z = %f, Y = %20.17e, Y_true = %20.17e\n", m, z, Y, Y_true[m + 3]);
                printf("Checking, abs tol = %e ... ", tol);
                if (abs(Y - Y_true[m + 3]) < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!\n");
                }
            }
        }

        printf("==============================================\n");
        printf("Test Mathieu modmc1\n");
        double mc1;
        double mc1d;
        N = 6;
        {
            double q = 0.01;
            double v = 4.0;
            double tol = 1e-12;
            double rtol = 1e-11;
            // Golden values obtained using Matlab.  I get values
            // for m = 1, 4, 7, 10, 13, 16.
            // I use widely spaced orders to verify adaption works OK.
            double mc1_true[6] = {-3.43416995264759106e-01, 3.97873130174148326e-01, 8.35567375774955712e-02,
                                  3.15168992451903213e-03,  4.36682920825898135e-05, 2.91722697679619490e-07};
            double mc1d_true[6] = {2.32628352485249795e-01, -1.37343635827296445e-01, 4.09262446440566696e-01,
                                   2.69766257109354118e-02, 5.19428246083634753e-04,  4.40526164453159552e-06};

            printf("q = %f\n", q);
            for (int i = 0; i < N; i++) {
                int m = 3 * i + 1;
                retcode = mathieu_modmc1(m, q, v, &mc1, &mc1d);
                printf(
                    "Iteration i = %d, order m = %d, v = %f, mc1 = %20.17e,      mc1d = %20.17e\n", i, m, v, mc1, mc1d
                );
                printf(
                    "                                        mc1_true = %20.17e, mc1d_true = %20.17e\n", mc1_true[i],
                    mc1d_true[i]
                );
                // For small 1 modmc1's output is very small.  To adequately test, I check
                // both the absolute err as well as the rel error.  I get one failure
                // which suggests some minute difference between Matlab and C imps.
                // The difference occurs because of catastrophic cancellation, but
                // I haven't found why Matlab and C are different yet.
                printf("Checking abs err, tol = %e ... ", tol);
                double diff = abs((mc1 - mc1_true[i]));
                double diffd = abs((mc1d - mc1d_true[i]));
                // printf("diff = %e, diffd = %e\n", diff, diffd);
                if (diff < tol && diffd < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  tol = %e, diff = %e, diffd = %e\n", tol, diff, diffd);
                }
                printf("Checking rel err, rtol = %e... ", rtol);
                double rdiff = abs((mc1 - mc1_true[i]) / mc1);
                double rdiffd = abs((mc1d - mc1d_true[i]) / mc1d);
                // printf("rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                if (rdiff < rtol && rdiffd < rtol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  rtol = %e, rdiff = %e, rdiffd = %e\n", rtol, rdiff, rdiffd);
                }
                printf("--------------------------------\n");
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 100;
            double v = 4.0;
            double tol = 1e-14;
            double rtol = 1e-13;
            // Golden values obtained using Matlab.  I get values
            // for m = 1, 4, 7, 10, 13, 16.
            // I use widely spaced orders to verify adaption works OK.
            double mc1_true[6] = {-3.41428322016207667e-02, 3.44948654167566600e-03,  3.35219048366175096e-02,
                                  -9.08021166319475555e-03, -3.24223338845559816e-02, 1.27646351014092524e-02};
            double mc1d_true[6] = {2.00691110743774284e-02,  1.85469526520107095e+01, -3.57255502109930578e+00,
                                   -1.79656659177568585e+01, 5.87447048444966935e+00, 1.72822185411075182e+01};

            printf("q = %f\n", q);
            for (int i = 0; i < N; i++) {
                int m = 3 * i + 1;
                retcode = mathieu_modmc1(m, q, v, &mc1, &mc1d);
                printf(
                    "Iteration i = %d, order m = %d, v = %f, mc1 = %20.17e,      mc1d = %20.17e\n", i, m, v, mc1, mc1d
                );
                printf(
                    "                                        mc1_true = %20.17e, mc1d_true = %20.17e\n", mc1_true[i],
                    mc1d_true[i]
                );
                printf("Checking abs err, tol = %e ... ", tol);
                double diff = abs((mc1 - mc1_true[i]));
                double diffd = abs((mc1d - mc1d_true[i]));
                // printf("diff = %e, diffd = %e\n", diff, diffd);
                if (diff < tol && diffd < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  tol = %e, diff = %e, diffd = %e\n", tol, diff, diffd);
                }
                printf("Checking rel err, rtol = %e... ", rtol);
                double rdiff = abs((mc1 - mc1_true[i]) / mc1);
                double rdiffd = abs((mc1d - mc1d_true[i]) / mc1d);
                // printf("rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                if (rdiff < rtol && rdiffd < rtol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  rtol = %e, rdiff = %e, rdiffd = %e\n", rtol, rdiff, rdiffd);
                }

                printf("--------------------------------\n");
            }
        }

        printf("==============================================\n");
        printf("Test Mathieu modms1\n");
        double ms1;
        double ms1d;
        N = 6;
        {
            double q = 0.01;
            double v = 4.0;
            double tol = 1e-6;
            double rtol = 1e-11; // I use different tol for relative err test.
            // Golden values obtained using Matlab.  I used
            // m = 1, 4, 7, 10, 13, 16
            double ms1_true[6] = {-3.43379253985569288e-01, 3.97873130174058565e-01, 8.35567375774849824e-02,
                                  3.15168992462108635e-03,  4.36682920825898135e-05, 2.91722697679619490e-07};
            double ms1d_true[6] = {2.29156633917084712e-01, -1.37343635825750598e-01, 4.09262446440550431e-01,
                                   2.69766257110715633e-02, 5.19428246083634753e-04,  4.40526164453159637e-06};

            printf("q = %f\n", q);
            for (int i = 0; i < N; i++) {
                int m = 3 * i + 1;
                retcode = mathieu_modms1(m, q, v, &ms1, &ms1d);
                printf(
                    "Iteration i = %d, order m = %d, v = %f, ms1 = %20.17e,      ms1d = %20.17e\n", i, m, v, ms1, ms1d
                );
                printf(
                    "                                        ms1_true = %20.17e, ms1d_true = %20.17e\n", ms1_true[i],
                    ms1d_true[i]
                );
                // For small 1 modmc1's output is very small.  To adequately test, I check
                // both the absolute err as well as the rel error.  I get one failure
                // which suggests some minute difference between Matlab and C imps.
                // The difference occurs because of catastrophic cancellation, but
                // I haven't found why Matlab and C are different yet.
                printf("Checking abs err, tol = %e ... ", tol);
                double diff = abs((ms1 - ms1_true[i]));
                double diffd = abs((ms1d - ms1d_true[i]));
                // printf("diff = %e, diffd = %e\n", diff, diffd);
                if (diff < tol && diffd < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  tol = %e, diff = %e, diffd = %e\n", tol, diff, diffd);
                }
                printf("Checking rel err, rtol = %e... ", rtol);
                double rdiff = abs((ms1 - ms1_true[i]) / ms1);
                double rdiffd = abs((ms1d - ms1d_true[i]) / ms1d);
                // printf("rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                if (rdiff < rtol && rdiffd < rtol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  rtol = %e, rdiff = %e, rdiffd = %e\n", rtol, rdiff, rdiffd);
                }
                printf("--------------------------------\n");
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 100;
            double v = 4.0;
            double tol = 1e-14;
            double rtol = 1e-14;
            // Golden values obtained using Matlab.  I used
            // m = 1, 4, 7, 10, 13, 16
            double ms1_true[6] = {-3.41201739333336501e-02, 2.33952744720816123e-03,  3.36955475281936756e-02,
                                  -8.29453993546177236e-03, -3.24933710575481122e-02, 1.27631677399401809e-02};
            double ms1d_true[6] = {-6.45251275907079758e-01, 1.85996271072809449e+01, -3.04135814724805620e+00,
                                   -1.80792391857584605e+01, 5.75576271594593436e+00, 1.72825425245382398e+01};

            printf("q = %f\n", q);
            for (int i = 0; i < N; i++) {
                int m = 3 * i + 1;
                retcode = mathieu_modms1(m, q, v, &ms1, &ms1d);
                printf(
                    "Iteration i = %d, order m = %d, v = %f, ms1 = %20.17e,      ms1d = %20.17e\n", i, m, v, ms1, ms1d
                );
                printf(
                    "                                        ms1_true = %20.17e, ms1d_true = %20.17e\n", ms1_true[i],
                    ms1d_true[i]
                );
                printf("Checking abs err, tol = %e ... ", tol);
                double diff = abs((ms1 - ms1_true[i]));
                double diffd = abs((ms1d - ms1d_true[i]));
                // printf("diff = %e, diffd = %e\n", diff, diffd);
                if (diff < tol && diffd < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  tol = %e, diff = %e, diffd = %e\n", tol, diff, diffd);
                }
                printf("Checking rel err, rtol = %e... ", rtol);
                double rdiff = abs((ms1 - ms1_true[i]) / ms1);
                double rdiffd = abs((ms1d - ms1d_true[i]) / ms1d);
                // printf("rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                if (rdiff < rtol && rdiffd < rtol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! rtol = %e, rdiff = %e, rdiffd = %e\n", rtol, rdiff, rdiffd);
                }
                printf("--------------------------------\n");
            }
        }

        printf("==============================================\n");
        printf("Test Mathieu modmc2\n");
        double mc2;
        double mc2d;
        N = 6;
        {
            double v = 4.0; // Test at v=4 since the value is extremely large for
                            // smaller v.
            double q = 0.01;
            double tol = 1e-9;
            double rtol = 1e-12; // I use different tol for relative err test.

            // Golden values obtained using Matlab.  I get values
            // for m = 1, 4, 7, 10, 13, 16.
            // I use widely spaced orders to verify adaption works OK.
            double mc2_true[6] = {-1.05364049460497881e-02, -5.43273917516489124e-02, -8.98079478272420850e-01,
                                  -1.21204986957549021e+01, -6.18466504094926336e+02, -7.25758230510118883e+04};
            double mc2d_true[6] = {-1.84664333620107901e+00, 1.61881073394890906e+00, 3.22020193455027437e+00,
                                   9.82487563907863120e+01,  7.22196427512374794e+03, 1.08632040247457405e+06};

            printf("q = %f\n", q);
            for (int i = 0; i < N; i++) {
                int m = 3 * i + 1;
                retcode = mathieu_modmc2(m, q, v, &mc2, &mc2d);
                printf(
                    "Iteration i = %d, order m = %d, v = %f, mc2 = %20.17e,      mc2d = %20.17e\n", i, m, v, mc2, mc2d
                );
                printf(
                    "                                        mc2_true = %20.17e, mc2d_true = %20.17e\n", mc2_true[i],
                    mc2d_true[i]
                );
                // This fcns output varies widely.  To adequately test, I check
                // both the absolute err as well as the rel error.
                printf("Checking abs err, tol = %e ... ", tol);
                double diff = abs((mc2 - mc2_true[i]));
                double diffd = abs((mc2d - mc2d_true[i]));
                // printf("diff = %e, diffd = %e\n", diff, diffd);
                if (diff < tol && diffd < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  tol = %e, diff = %e, diffd = %e\n", tol, diff, diffd);
                }
                printf("Checking rel err, rtol = %e... ", rtol);
                double rdiff = abs((mc2 - mc2_true[i]) / mc2);
                double rdiffd = abs((mc2d - mc2d_true[i]) / mc2d);
                // printf("rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                if (rdiff < rtol && rdiffd < rtol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! rtol = %e,  rdiff = %e, rdiffd = %e\n", rtol, rdiff, rdiffd);
                }
                printf("--------------------------------\n");
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 100;
            double v = 4.0; // Test at v=4 since the value is extremely large for smaller v.
            double tol = 1e-14;
            double rtol = 1e-12;
            // Golden values obtained using Matlab.  I get values
            // for m = 1, 4, 7, 10, 13, 16.
            // I use widely spaced orders to verify adaption works OK.
            double mc2_true[6] = {-5.50400810104661117e-06, -3.39713077181867523e-02, 6.51339805952160303e-03,
                                  3.29219862870713670e-02,  -1.07334892119392861e-02, -3.16798833654585860e-02};
            double mc2d_true[6] = {-1.86457777769476074e+01, 1.90043837315601305e+00,  1.82970001977819763e+01,
                                   -4.97283180832788130e+00, -1.76904663568990372e+01, 6.98187640938441945e+00};

            printf("q = %f\n", q);
            for (int i = 0; i < N; i++) {
                int m = 3 * i + 1;
                retcode = mathieu_modmc2(m, q, v, &mc2, &mc2d);
                printf(
                    "Iteration i = %d, order m = %d, v = %f, mc2 = %20.17e,      mc2d = %20.17e\n", i, m, v, mc2, mc2d
                );
                printf(
                    "                                        mc2_true = %20.17e, mc2d_true = %20.17e\n", mc2_true[i],
                    mc2d_true[i]
                );
                printf("Checking abs err, tol = %e ... ", tol);
                double diff = abs((mc2 - mc2_true[i]));
                double diffd = abs((mc2d - mc2d_true[i]));
                // printf("diff = %e, diffd = %e\n", diff, diffd);
                if (diff < tol && diffd < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  diff = %e, diffd = %e\n", diff, diffd);
                }
                printf("Checking rel err, rtol = %e... ", rtol);
                double rdiff = abs((mc2 - mc2_true[i]) / mc2);
                double rdiffd = abs((mc2d - mc2d_true[i]) / mc2d);
                // printf("rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                if (rdiff < rtol && rdiffd < rtol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                }
                printf("--------------------------------\n");
            }
        }

        printf("==============================================\n");
        printf("Test Mathieu modms2\n");
        double ms2;
        double ms2d;
        N = 6;
        {
            double q = 0.01;
            double v = 4.0; // Test at v=4 since the value is extremely large for smaller v.
            double tol = 1e-2;
            double rtol = 1e-8;
            // Golden values obtained using Matlab.  I used
            // m = 1, 4, 7, 10, 13, 16
            double ms2_true[6] = {
                -9.91337437572497975e-03, -5.43273917519763866e-02, -8.98079478272391984e-01,
                -1.21204986957547884e+01, -6.18466504094926336e+02, -7.25758230510119029e+04,
            };
            double ms2d_true[6] = {
                -1.84736861502861638e+00, 1.61881073394918618e+00, 3.22020193455016246e+00,
                9.82487563907851893e+01,  7.22196427512374976e+03, 1.08632040247457428e+06,
            };

            printf("q = %f\n", q);
            for (int i = 0; i < N; i++) {
                int m = 3 * i + 1;
                retcode = mathieu_modms2(m, q, v, &ms2, &ms2d);
                printf(
                    "Iteration i = %d, order m = %d, v = %f, ms2 = %20.17e,      ms2d = %20.17e\n", i, m, v, ms2, ms2d
                );
                printf(
                    "                                        ms2_true = %20.17e, ms2d_true = %20.17e\n", ms2_true[i],
                    ms2d_true[i]
                );
                // To adequately test, I check
                // both the absolute err as well as the rel error.
                printf("Checking abs err, tol = %e ... ", tol);
                double diff = abs((ms2 - ms2_true[i]));
                double diffd = abs((ms2d - ms2d_true[i]));
                // printf("diff = %e, diffd = %e\n", diff, diffd);
                if (diff < tol && diffd < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  tol = %e, diff = %e, diffd = %e\n", tol, diff, diffd);
                }
                printf("Checking rel err, rtol = %e... ", rtol);
                double rdiff = abs((ms2 - ms2_true[i]) / ms2);
                double rdiffd = abs((ms2d - ms2d_true[i]) / ms2d);
                // printf("rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                if (rdiff < rtol && rdiffd < rtol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  rtol = %e, rdiff = %e, rdiffd = %e\n", rtol, rdiff, rdiffd);
                }
                printf("--------------------------------\n");
            }
        }
        printf("----------------------------------------------\n");
        {
            double q = 100;
            double v = 4.0; // Test at v=4 since the value is extremely large for smaller v.
            double tol = 1e-14;
            double rtol = 1e-14;
            // Golden values obtained using Matlab.  I used
            // m = 1, 4, 7, 10, 13, 16
            double ms2_true[6] = {1.21267949732755992e-03, -3.40647275296700053e-02, 5.53991666556525025e-03,
                                  3.31278916985982691e-02, -1.05157990917375058e-02, -3.16804730008901705e-02};
            double ms2d_true[6] = {-1.86352300142727003e+01, 1.29451048253596168e+00,  1.83932580755102961e+01,
                                   -4.54427791318153673e+00, -1.77295648089783491e+01, 6.98107651232898974e+00};

            printf("q = %f\n", q);
            for (int i = 0; i < N; i++) {
                int m = 3 * i + 1;
                retcode = mathieu_modms2(m, q, v, &ms2, &ms2d);
                printf(
                    "Iteration i = %d, order m = %d, v = %f, ms2 = %20.17e,      ms2d = %20.17e\n", i, m, v, ms2, ms2d
                );
                printf(
                    "                                        ms2_true = %20.17e, ms2d_true = %20.17e\n", ms2_true[i],
                    ms2d_true[i]
                );
                printf("Checking abs err, tol = %e ... ", tol);
                double diff = abs((ms2 - ms2_true[i]));
                double diffd = abs((ms2d - ms2d_true[i]));
                // printf("diff = %e, diffd = %e\n", diff, diffd);
                if (diff < tol && diffd < tol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed!  tol = %e, diff = %e, diffd = %e\n", tol, diff, diffd);
                }
                printf("Checking rel err, rtol = %e... ", rtol);
                double rdiff = abs((ms2 - ms2_true[i]) / ms2);
                double rdiffd = abs((ms2d - ms2d_true[i]) / ms2d);
                // printf("rdiff = %e, rdiffd = %e\n", rdiff, rdiffd);
                if (rdiff < rtol && rdiffd < rtol) {
                    pass++;
                    printf("Passed!\n");
                } else {
                    fail++;
                    printf("Failed! rtol = %e, rdiff = %e, rdiffd = %e\n", rtol, rdiff, rdiffd);
                }
                printf("--------------------------------\n");
            }
        }

        printf("===========================================\n");
        printf("At end of run, pass = %d, fail = %d\n", pass, fail);

        return 0;
    }

} // namespace mathieu
} // namespace xsf
