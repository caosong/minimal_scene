/* kruppa.c */
/* Decompose an F-matrix using the Krupa equations */

#include <float.h>
#include <stdio.h>

#include "matrix.h"

#include "kruppa.h"

int decompose_fmatrix_kruppa(double *F, double *focal)
{
    double f_min = 0.05;
    double f_max = 10.0;
    double f_step = 0.01;

    double v1_sq[3], v2_sq[3], v12[3], u1_sq[3], u2_sq[3], u12[3];
    double f, f_best = 0.0;
    double min_error = DBL_MAX;
    
    /* Run SVD on F */
    double U[9], VT[9], S[3];
    dgesvd_driver(3, 3, F, U, S, VT);

    v1_sq[0] = VT[0] * VT[0];
    v1_sq[1] = VT[1] * VT[1];
    v1_sq[2] = VT[2] * VT[2];

    v2_sq[0] = VT[3] * VT[3];
    v2_sq[1] = VT[4] * VT[4];
    v2_sq[2] = VT[5] * VT[5];

    v12[0] = VT[0] * VT[3];
    v12[1] = VT[1] * VT[4];
    v12[2] = VT[2] * VT[5];

    u1_sq[0] = U[0] * U[0];
    u1_sq[1] = U[3] * U[3];
    u1_sq[2] = U[6] * U[6];

    u2_sq[0] = U[3] * U[3];
    u2_sq[1] = U[4] * U[4];
    u2_sq[2] = U[5] * U[5];

    u12[0] = U[0] * U[3];
    u12[1] = U[1] * U[4];
    u12[2] = U[2] * U[5];

    for (f = f_min; f <= f_max; f += f_step) {
        /* Compute error */
        double fvec[3] = { f, f, 1.0 };
        double n1, n2, n3, d1, d2, d3, e1, e2, e3, error;

        matrix_product(1, 3, 3, 1, v2_sq, fvec, &n1);
        matrix_product(1, 3, 3, 1, v12, fvec, &n2);
        matrix_product(1, 3, 3, 1, v1_sq, fvec, &n3);

        matrix_product(1, 3, 3, 1, u1_sq, fvec, &d1);
        matrix_product(1, 3, 3, 1, u12, fvec, &d2);
        matrix_product(1, 3, 3, 1, u2_sq, fvec, &d3);

        n2 = -n2;
        d1 = S[0] * S[0] * d1;
        d2 = S[0] * S[1] * d2;
        d3 = S[1] * S[1] * d3;

        e1 = n1 / d1;
        e2 = n2 / d2;
        e3 = n3 / d3;

        error = (e1 - e2) * (e1 - e2) + (e1 - e3) * (e1 - e3) + 
            (e2 - e3) * (e2 - e3);
        
        if (error < min_error) {
            min_error = error;
            f_best = f;
        }
    }

    *focal = f_best;
    printf("[decompose_fmatrix_kruppa] "
           "Best focal length: %0.3f, error: %0.3f\n", f_best, min_error);

    return 0;
}
