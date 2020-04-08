//
// Created by vover on 3/11/20.
//

#include <stdlib.h>
#include <stdio.h>
#include "harmonic_help.h"
#include "harmonic.h"
#include "residuals.h"
#include "laspack/qmatrix.h"
#include "laspack/rtc.h"
#include "case.h"
#include "differential_operators.h"

#define N_X 30
#define N_Y 30
#define L_X 2*M_PI
#define L_Y 2*M_PI

typedef struct {
    double C_mn;
    double D_mn;
    double P_mn;
} coeff;

void arrays_init (double **u_out, double **u_in, double **u_curr, double **u_new,
                  double **v_out, double **v_in, double **v_curr, double **v_new,
                  double **p_out, double **p_in, double **p_curr, double **p_new, int N_x, int N_y);
void arrays_free (double **u_out, double **u_in, double **u_curr, double **u_new,
                  double **v_out, double **v_in, double **v_curr, double **v_new,
                  double **p_out, double **p_in, double **p_curr, double **p_new);
void array_zero (double *a, int size);
void array_to_array (double *from, double *to, int size);

void solve_first () {

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

    eigenvalue_mn (&eval, p_g.mu, m, n);
    eigenvector_mn (&evec, p_g.mu, m, n);

    fill_with_vector(G, V1, V2, p_s, m, n, evec.v3);

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
        coef.C_mn = coefficient_Cmn(V1, p_s, m, n);
        coef.D_mn = coefficient_Dmn(V2, p_s, m, n);
        coef.P_mn = coefficient_Pmn(G, p_s, m, n);

        printf("\n %lf \t %lf \t %lf  \n", coef.C_mn, coef.D_mn, coef.P_mn);
        printf("\n %lf \t %lf \t %lf  \n", coef.C_mn/coeff_f.C_mn, coef.D_mn/coeff_f.D_mn, coef.P_mn/coeff_f.P_mn);
    }


    return;
}

void solve_second() {
    double *u_out, *u_in, *u_curr, *u_new;
    double *v_out, *v_in, *v_curr, *v_new;
    double *p_out, *p_in, *p_curr, *p_new;
    int N_x = N_X, N_y = N_Y;
    double L_x = L_X, L_y = L_Y;
    double h_x = L_X / N_X, h_y = L_Y / N_Y, tay = 0.01;
    int u_size = (N_X + 1) * N_Y;
    int v_size = (N_Y + 1) * N_X;
    int p_size = N_X * N_Y;
    double p_ro = 1;
    double mu = 0.1;
    double c_ro_0 = 1;
    double ro__0 = 0.01;
    int m=1,n=1, j, k;
    eigenvalue eval;
    eigenvector evec;
    coeff coeff_f, coef;

    eigenvalue_mn (&eval, mu, m, n);
    eigenvector_mn (&evec, mu, m, n);


    printf("vector 1: (%lf, \t %lf, \t %lf)\n"
           "vector 2: (%lf, \t %lf, \t %lf) \n"
           "vector 3: (%lf, \t %lf, \t %lf) \n",
           evec.v1[0], evec.v1[1], evec.v1[2],
           evec.v2[0], evec.v2[1], evec.v2[2],
           evec.v3[0], evec.v3[1], evec.v3[2] );

    arrays_init (&u_out, &u_in, &u_curr, &u_new,
                      &v_out, &v_in, &v_curr, &v_new,
                      &p_out, &p_in, &p_curr, &p_new, N_x, N_y);

    fill_with_vector_new (u_curr, v_curr, p_curr, N_x, N_y, h_x, h_y, m, n, evec.v2);


    coeff_f.C_mn = coefficient_Cmn_new(u_curr, u_size, N_x, N_y, h_x, h_y, m, n);
    coeff_f.D_mn = coefficient_Dmn_new(v_curr, v_size, N_x, N_y, h_x, h_y, m, n);
    coeff_f.P_mn = coefficient_Pmn_new(p_curr, p_size, N_x, N_y, h_x, h_y, m, n);


    printf("\n %lf \t %lf \t %lf \n",
           coeff_f.C_mn, coeff_f.D_mn, coeff_f.P_mn);

    for (int i = 0; i < 3; ++i){

        array_to_array(u_curr, u_in, u_size);
        array_to_array(v_curr, v_in, v_size);
        array_to_array(p_curr, p_in, p_size);

        udxdx(u_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
        vdydy(v_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
        for (j = 0; j < u_size; ++j)
            u_new[j] = mu * tay * 4 * u_out[j]/ (ro__0 * 3);
        for (k = 0; k < v_size; ++k)
            v_new[k] = mu * tay * 4 * v_out[k]/ (ro__0 * 3);

        udydy(u_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
        vdxdx(v_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
        for (j = 0; j < u_size; ++j)
            u_new[j] += mu * tay * u_out[j]/ (ro__0);
        for (k = 0; k < v_size; ++k)
            v_new[k] += mu * tay * v_out[k]/ (ro__0);

        udxdy(u_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
        vdxdy(v_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
        for (j = 0; j < u_size; ++j)
            u_new[j] += mu * tay * u_out[j]/ (ro__0 * 3);
        for (k = 0; k < v_size; ++k)
            v_new[k] += mu * tay * v_out[k]/ (ro__0 * 3);

        pdx(u_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
        pdy(v_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
        for (j = 0; j < u_size; ++j)
            u_new[j] -= c_ro_0 * tay * u_out[j]/ (ro__0);
        for (k = 0; k < v_size; ++k)
            v_new[k] -= c_ro_0 * tay * v_out[k]/ (ro__0);

        for (j = 0; j < u_size; ++j)
            u_new[j] += u_in[j];
        for (k = 0; k < v_size; ++k)
            v_new[k] += v_in[k];

        udx(p_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
        for (j = 0; j < p_size; ++j)
            p_new[j] = -p_out[j] * ro__0 * tay;
        vdy(p_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
        for (j = 0; j < p_size; ++j)
            p_new[j] += -p_out[j] * ro__0 * tay;

        for (j = 0; j < p_size; ++j)
            p_new[j] += p_in[j];

        coef.C_mn = coefficient_Cmn_new(u_new, u_size, N_x, N_y, h_x, h_y, m, n);
        coef.D_mn = coefficient_Dmn_new(v_new, v_size, N_x, N_y, h_x, h_y, m, n);
        coef.P_mn = coefficient_Pmn_new(p_new, p_size, N_x, N_y, h_x, h_y, m, n);

        printf("\n %lf \t %lf \t %lf  \t %lf \n",
               coef.C_mn, coef.D_mn, coef.P_mn, residial_Ch_(u_new,v_new,p_new,u_size, v_size,p_size));
        printf("\n %lf \t %lf \t %lf  \n", coef.C_mn/coeff_f.C_mn, coef.D_mn/coeff_f.D_mn, coef.P_mn/coeff_f.P_mn);

        array_to_array(u_new, u_curr, u_size);
        array_to_array(v_new, v_curr, v_size);
        array_to_array(p_new, p_curr, p_size);
    }


    arrays_free (&u_out, &u_in, &u_curr, &u_new,
                      &v_out, &v_in, &v_curr, &v_new,
                      &p_out, &p_in, &p_curr, &p_new);


    return;
}

void arrays_init (double **u_out, double **u_in, double **u_curr, double **u_new,
                  double **v_out, double **v_in, double **v_curr, double **v_new,
                  double **p_out, double **p_in, double **p_curr, double **p_new, int N_x, int N_y){

    *u_out = (double *) malloc (sizeof(double) * (N_x + 1) * N_y);
    *u_in  = (double *) malloc (sizeof(double) * (N_x + 1) * N_y);
    *u_curr= (double *) malloc (sizeof(double) * (N_x + 1) * N_y);
    *u_new = (double *) malloc (sizeof(double) * (N_x + 1) * N_y);

    *v_out = (double *) malloc (sizeof(double) * (N_y + 1) * N_x);
    *v_in  = (double *) malloc (sizeof(double) * (N_y + 1) * N_x);
    *v_curr= (double *) malloc (sizeof(double) * (N_y + 1) * N_x);
    *v_new = (double *) malloc (sizeof(double) * (N_y + 1) * N_x);

    *p_out = (double *) malloc (sizeof(double) * N_x * N_y);
    *p_in  = (double *) malloc (sizeof(double) * N_x * N_y);
    *p_curr= (double *) malloc (sizeof(double) * N_x * N_y);
    *p_new = (double *) malloc (sizeof(double) * N_x * N_y);


    return;
}
void arrays_free (double **u_out, double **u_in, double **u_curr, double **u_new,
                  double **v_out, double **v_in, double **v_curr, double **v_new,
                  double **p_out, double **p_in, double **p_curr, double **p_new){

    free(*u_out); free(*u_in); free(*u_curr); free(*u_new);
    free(*v_out); free(*v_in); free(*v_curr); free(*v_new);
    free(*p_out); free(*p_in); free(*p_curr); free(*p_new);

    return;
}

void array_zero (double *a, int size) {

    for (int i = 0; i < size; i++)
        a[i] = 0;

    return;
}

void array_to_array (double *from, double *to, int size) {

    for (int i = 0; i < size; ++i)
        to[i] = from[i];

    return;
}