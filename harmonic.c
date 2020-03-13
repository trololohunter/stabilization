//
// Created by vover on 3/11/20.
//

#include <stdlib.h>
#include <stdio.h>
#include "harmonic_help.h"
#include "residuals.h"
#include "laspack/qmatrix.h"
#include "laspack/rtc.h"
#include "case.h"

typedef struct {
    double C_mn;
    double D_mn;
    double P_mn;
} coeff;

void solve () {

    P_gas p_g;
    P_she p_s;
    double *G;
    double *V1;
    double *V2;
    int *st;
    int it_sp=1, it_t=1, k;
    Norm_Step n_s, n_first;
    QMatrix_L A;
    Vector b, x;
    T_const t_c;
    MUM_const m_c;
    MM_step m_s;
    double GG = 0;
    double tmp;
    size_t mm;
    int m=1,n=1,mmm;
    eigenvalue eval;
    eigenvector evec;
    coeff coeff_f, coef;


    param_dif(&p_g);
    param_she_step(&p_s, p_g, it_t, it_sp);

    st = (int*) malloc(p_s.Dim * sizeof(int));
    G = (double*) malloc(p_s.Dim * sizeof(double));
    V1 = (double*) malloc(p_s.Dim * sizeof(double));
    V2 = (double*) malloc(p_s.Dim * sizeof(double));

    Setka(st, &p_s);

    eigenvalue_mn (&eval, p_g, m, n);
    eigenvector_mn (&evec, p_g, m, n);

    fill_with_vector(G, V1, V2, p_s, m, n, evec.v1);

    param_t_const(&t_c, p_s, p_g);
    SetRTCAccuracy(1e-8);

    coeff_f.C_mn = coefficient_Cmn(V1, p_s, m, n);
    coeff_f.D_mn = coefficient_Dmn(V2, p_s, m, n);
    coeff_f.P_mn = coefficient_Pmn(G, p_s, m, n);

    printf("\n %lf \t %lf \t %lf  \n", coeff_f.C_mn, coeff_f.D_mn, coeff_f.P_mn);

    for (k = 1; k < p_s.N + 1; ++k)
    {
        mm = 1;

        Q_Constr(&A, "A", (size_t) 3 * p_s.Dim, False, Rowws, Normal, True);
        V_Constr(&b, "b", (size_t) 3 * p_s.Dim, Normal, True);
        V_Constr(&x, "x", (size_t) 3 * p_s.Dim, Normal, True);

        GG = -1000000;
        for (mmm = 0; mmm < p_s.Dim; ++mmm)
        {
//            printf("%lf \n", GG);
            if (exp(-G[mmm]) > GG) GG = exp(-G[mmm]);
        }
        //printf("%lf \n", GG);
        param_MUM_const(&m_c, p_s, GG, p_g);

        for (size_t i = 0; i < p_s.Dim; ++i)
        {
            V_SetCmp(&x, 3 * i + 1, G[i]);
            V_SetCmp(&x, 3 * i + 2, V1[i]);
            V_SetCmp(&x, 3 * i + 3, V2[i]);
        }


        for (mmm = 0; mmm < p_s.Dim; ++mmm)
        {
            param_MM_step(&m_s,mm,p_s.M_x+1, mmm);
            switch (st[mmm])
            {
                case 0:
                    mm = case0(&A, &b, t_c, m_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, p_g.p_ro, mm);
                    break;
                case 1:
                    mm = case1(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm);
                    break;
                case 2:
                    mm = case2(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm);
                    break;
                case 3:
                    mm = case3(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm, p_g.omega);
                    break;
                case 4:
                    mm = case4(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm);
                    break;
                case 5:
                    mm = case1(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm);
                    break;
                case 6:
                    mm = case2(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm);
                    break;
                case 7:
                    mm = case1(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm);
                    break;
                case 8:
                    mm = case2(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm);
                    break;
                case 9:
                    if (WALL) mm = case9(&A, &b, t_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, mm);
                    else mm = case0(&A, &b, t_c, m_c, m_s, k, V1, V2, G, mmm, p_s, p_g.mu, p_g.p_ro, mm);
                    break;
                case 10:
                    mm = case10(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_g.mu, mm);
                    break;
                default:
                {
                    printf("FATAL ERROR: unknown type of node");
                    exit(1);
                }

            }
            // printf("%d \n", m);
            ++mm;
        }

        CGSIter(&A, &x, &b, 2000, SSORPrecond, 1);

        for (size_t i = 0; i < p_s.Dim; ++i)
        {
            G[i] = V_GetCmp(&x, 3 * i + 1);
            V1[i] = V_GetCmp(&x, 3 * i + 1 + 1);
            V2[i] = V_GetCmp(&x, 3 * i + 2 + 1);
        }

        for (mmm = 0; mmm < p_s.Dim; ++mmm)
        {
            if (fabs(G[mmm])  < 1e-14 ) G[mmm] = 0;
            if (fabs(V1[mmm]) < 1e-14 ) V1[mmm] = 0;
            if (fabs(V2[mmm]) < 1e-14 ) V2[mmm] = 0;
        }
        //usleep(10);


        Q_Destr(&A);
        V_Destr(&b);
        V_Destr(&x);

        residual_Ch(V1, V2, G, p_s, &n_s, u1, u2, g);
        //fprintf(fp, "%lf \t %lf \n", k*p_s.tau, sqrt(n_s.V1norm * n_s.V1norm + n_s.V2norm * n_s.V2norm));
        //if (SMOOTH_SOLUTION != 1)
    }

    coef.C_mn = coefficient_Cmn(V1, p_s, m, n);
    coef.D_mn = coefficient_Dmn(V2, p_s, m, n);
    coef.P_mn = coefficient_Pmn(G, p_s, m, n);

    printf("\n %lf \t %lf \t %lf  \n", coef.C_mn, coef.D_mn, coef.P_mn);
    printf("\n %lf \t %lf \t %lf  \n", coef.C_mn/coeff_f.C_mn, coef.D_mn/coeff_f.D_mn, coef.P_mn/coeff_f.P_mn);


    return;
}