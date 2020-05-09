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
#define N_X_ 50
#define N_Y_ 50

typedef struct {
    double C_mn;
    double D_mn;
    double P_mn;
} coeff;

void arrays_init_for_second (double **u_out, double **u_in, double **u_curr, double **u_new,
                  double **v_out, double **v_in, double **v_curr, double **v_new,
                  double **p_out, double **p_in, double **p_curr, double **p_new, int N_x, int N_y);
void arrays_free_for_second (double **u_out, double **u_in, double **u_curr, double **u_new,
                  double **v_out, double **v_in, double **v_curr, double **v_new,
                  double **p_out, double **p_in, double **p_curr, double **p_new);
void arrays_init_for_third (double **u_out, double **u_in, double **u_curr, double **u_new,
                             double **v_out, double **v_in, double **v_curr, double **v_new,
                             double **p_out, double **p_in, double **p_curr, double **p_new,
                            double **P_in, double **P_out, double **P,
                            double **C_in, double **C_out, double **C,
                            double **D_in, double **D_out, double **D,
                            int N_x, int N_y);
void arrays_free_for_third (double **u_out, double **u_in, double **u_curr, double **u_new,
                             double **v_out, double **v_in, double **v_curr, double **v_new,
                             double **p_out, double **p_in, double **p_curr, double **p_new,
                            double **P_in, double **P_out, double **P,
                            double **C_in, double **C_out, double **C,
                            double **D_in, double **D_out, double **D);
void array_zero (double *a, int size);
void array_to_array (double *from, double *to, int size);
void print_matr (double *a, int N, int M);
void print_evals (double *eval, int *ind_m, int *ind_n, int *ind_l);
void sort_evals (double *eval, int *ind_m, int *ind_n, int *ind_l);

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
    arrays_init_for_second(&u_out, &u_in, &u_curr, &u_new,
                           &v_out, &v_in, &v_curr, &v_new,
                           &p_out, &p_in, &p_curr, &p_new, N_x, N_y);
for (m = 1; m < 10; m++)
    for (n = 1; n < 10; n++) {


        fill_with_vector_new(u_curr, v_curr, p_curr, N_x, N_y, h_x, h_y, 1, 1, evec.v2);

        //print_matr(u_curr, N_y, N_x + 1);
        //print_matr(v_curr, N_y + 1, N_x);
        //print_matr(p_curr, N_y, N_x);

        coeff_f.C_mn = coefficient_Cmn_new(u_curr, u_size, N_x, N_y, h_x, h_y, m, n);
        coeff_f.D_mn = coefficient_Dmn_new(v_curr, v_size, N_x, N_y, h_x, h_y, m, n);
        coeff_f.P_mn = coefficient_Pmn_new(p_curr, p_size, N_x, N_y, h_x, h_y, m, n);


        printf("\n C_%d%d = %lf \t D_%d%d = %lf \t P_%d%d = %lf \n",
               m,n, coeff_f.C_mn, m,n,coeff_f.D_mn, m,n,coeff_f.P_mn);

        for (int i = 0; i < 3; ++i) {

            array_to_array(u_curr, u_in, u_size);
            array_to_array(v_curr, v_in, v_size);
            array_to_array(p_curr, p_in, p_size);

            udxdx(u_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
            vdydy(v_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
            for (j = 0; j < u_size; ++j)
                u_new[j] = mu * tay * 4 * u_out[j] / (ro__0 * 3);
            for (k = 0; k < v_size; ++k)
                v_new[k] = mu * tay * 4 * v_out[k] / (ro__0 * 3);

            udydy(u_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
            vdxdx(v_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
            for (j = 0; j < u_size; ++j)
                u_new[j] += mu * tay * u_out[j] / (ro__0);
            for (k = 0; k < v_size; ++k)
                v_new[k] += mu * tay * v_out[k] / (ro__0);

            udxdy(u_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
            vdxdy(v_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
            for (j = 0; j < u_size; ++j)
                u_new[j] += mu * tay * u_out[j] / (ro__0 * 3);
            for (k = 0; k < v_size; ++k)
                v_new[k] += mu * tay * v_out[k] / (ro__0 * 3);

            pdx(u_out, u_in, N_x, N_y, h_x, h_y, L_x, L_y);
            pdy(v_out, v_in, N_x, N_y, h_x, h_y, L_x, L_y);
            for (j = 0; j < u_size; ++j)
                u_new[j] -= c_ro_0 * tay * u_out[j] / (ro__0);
            for (k = 0; k < v_size; ++k)
                v_new[k] -= c_ro_0 * tay * v_out[k] / (ro__0);

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

            //print_matr(u_new, N_y, N_x + 1);
            //print_matr(v_new, N_y + 1, N_x);
            //print_matr(p_new, N_y, N_x);

            coef.C_mn = coefficient_Cmn_new(u_new, u_size, N_x, N_y, h_x, h_y, m, n);
            coef.D_mn = coefficient_Dmn_new(v_new, v_size, N_x, N_y, h_x, h_y, m, n);
            coef.P_mn = coefficient_Pmn_new(p_new, p_size, N_x, N_y, h_x, h_y, m, n);

            printf("\n C_%d%d = %lf \t D_%d%d = %lf \t P_%d%d = %lf \t %lf\n",
                   m,n, coeff_f.C_mn, m,n,coeff_f.D_mn, m,n,coeff_f.P_mn, residial_Ch_(u_new, v_new, p_new, u_size, v_size, p_size));
            //printf("\n %lf \t %lf \t %lf  \t %lf \n",
            //       coef.C_mn, coef.D_mn, coef.P_mn, residial_Ch_(u_new, v_new, p_new, u_size, v_size, p_size));
            //printf("\n %lf \t %lf \t %lf  \n", coef.C_mn / coeff_f.C_mn, coef.D_mn / coeff_f.D_mn,
            //       coef.P_mn / coeff_f.P_mn);

            array_to_array(u_new, u_curr, u_size);
            array_to_array(v_new, v_curr, v_size);
            array_to_array(p_new, p_curr, p_size);
        }

    }
    arrays_free_for_second (&u_out, &u_in, &u_curr, &u_new,
                      &v_out, &v_in, &v_curr, &v_new,
                      &p_out, &p_in, &p_curr, &p_new);


    return;
}


void solve_third() {
    double *u_out, *u_in, *u_curr, *u_new;
    double *v_out, *v_in, *v_curr, *v_new;
    double *p_out, *p_in, *p_curr, *p_new;
    double *P_in, *P_out, *P;
    double *C_in, *C_out, *C;
    double *D_in, *D_out, *D;
    int N_x = N_X, N_y = N_Y;
    double L_x = L_X, L_y = L_Y;
    double h_x = L_X / N_X, h_y = L_Y / N_Y, tau = 0.01;
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
    arrays_init_for_third(&u_out, &u_in, &u_curr, &u_new,
                           &v_out, &v_in, &v_curr, &v_new,
                           &p_out, &p_in, &p_curr, &p_new,
                           &P_in, &P_out, &P,
                           &C_in, &C_out, &C,
                           &D_in, &D_out, &D,
                          N_x, N_y);
    for (m = 1; m < 2; m++)
        for (n = 1; n < 2; n++) {


            fill_with_vector_new(u_curr, v_curr, p_curr, N_x, N_y, h_x, h_y, 1, 1, evec.v2);

            //print_matr(u_curr, N_y, N_x + 1);
            //print_matr(v_curr, N_y + 1, N_x);
            //print_matr(p_curr, N_y, N_x);



            coeff_f.C_mn = coefficient_Cmn_new(u_curr, u_size, N_x, N_y, h_x, h_y, m, n);
            coeff_f.D_mn = coefficient_Dmn_new(v_curr, v_size, N_x, N_y, h_x, h_y, m, n);
            coeff_f.P_mn = coefficient_Pmn_new(p_curr, p_size, N_x, N_y, h_x, h_y, m, n);


            printf("\n C_%d%d = %lf \t D_%d%d = %lf \t P_%d%d = %lf \n",
                   m,n, coeff_f.C_mn, m,n,coeff_f.D_mn, m,n,coeff_f.P_mn);

            coef_from_function(u_curr, v_curr, p_curr, C, D, P, h_x, h_y, N_x, N_y);

            for (int i = 0; i < 3; ++i) {

                //array_to_array(u_curr, u_in, u_size);
                //array_to_array(v_curr, v_in, v_size);
                //array_to_array(p_curr, p_in, p_size);



                print_matr(C, (N_y-1), (N_x-1));
                print_matr(D, (N_y-1), (N_x-1));
                print_matr(P, (N_y-1), (N_x-1));

                array_to_array(P, P_in, N_x * N_y);
                array_to_array(C, C_in, N_x * N_y);
                array_to_array(D, D_in, N_x * N_y);

                count_next_coef(C_in, D_in, P_in, C_out, D_out, P_out, N_x, N_y, h_x, h_y, ro__0, c_ro_0, tau, mu);

            /*    if (i == 0) {
                    for (int ii = 0; ii < N_x; ++ii) {
                        //C[ii] = 0;
                        C[ii*N_x] = 0;
                    }
                    for (int jj = 0; jj < N_y; ++jj) {
                        //D[jj*N_x] = 0;
                        D[jj] = 0;
                    }
                }*/

                funcion_from_coef(C_out, D_out, P_out, u_new, v_new, p_new, h_x, h_y, N_x, N_y, mu, i*tau);

                coef.C_mn = coefficient_Cmn_new(u_new, u_size, N_x, N_y, h_x, h_y, m, n);
                coef.D_mn = coefficient_Dmn_new(v_new, v_size, N_x, N_y, h_x, h_y, m, n);
                coef.P_mn = coefficient_Pmn_new(p_new, p_size, N_x, N_y, h_x, h_y, m, n);

                printf("\n C_%d%d = %lf \t D_%d%d = %lf \t P_%d%d = %lf \t %e\n",
                       m,n, coeff_f.C_mn, m,n,coeff_f.D_mn, m,n,coeff_f.P_mn, residial_Ch_(u_new, v_new, p_new, u_size, v_size, p_size));

                print_matr(u_curr, N_y, N_x + 1);
                print_matr(v_curr, N_y + 1, N_x);
                print_matr(p_curr, N_y, N_x);

                array_to_array(P_out, P, N_x * N_y);
                array_to_array(C_out, C, N_x * N_y);
                array_to_array(D_out, D, N_x * N_y);

                array_to_array(u_new, u_curr, u_size);
                array_to_array(v_new, v_curr, v_size);
                array_to_array(p_new, p_curr, p_size);
            }

        }
    arrays_free_for_third (&u_out, &u_in, &u_curr, &u_new,
                            &v_out, &v_in, &v_curr, &v_new,
                            &p_out, &p_in, &p_curr, &p_new,
                           &P_in, &P_out, &P,
                           &C_in, &C_out, &C,
                           &D_in, &D_out, &D);


    return;
}

void qwerty () {
    eigenvalue eval;
    double *values = (double*) malloc(sizeof(double) * (N_X_-1)*(N_Y_-1)*3);
    int *indexes_m = (int*) malloc (sizeof(int) * (N_X_-1)*(N_Y_-1)*3);
    int *indexes_n = (int*) malloc (sizeof(int) * (N_X_-1)*(N_Y_-1)*3);
    int *indexes_l = (int*) malloc (sizeof(int) * (N_X_-1)*(N_Y_-1)*3);

    for (int n = 1; n < N_Y_; ++n)
        for (int m = 1; m < N_X_; ++m)
        {
            eigenvalue_mn(&eval, 0.1, m, n);
            printf ("%d \t %d \n", m, n);
            printf ("%lf \t %lf \t %lf \n", eval.lambda_1, eval.lambda_2, eval.lambda_3);
            values [((n-1)*(N_X_-1) + m-1)*3] = eval.lambda_1;
            indexes_m[((n-1)*(N_X_-1) + m-1)*3] = m;
            indexes_n[((n-1)*(N_X_-1) + m-1)*3] = n;
            indexes_l[((n-1)*(N_X_-1) + m-1)*3] = 1;
            values [((n-1)*(N_X_-1) + m-1)*3+1] = eval.lambda_2;
            indexes_m[((n-1)*(N_X_-1) + m-1)*3+1] = m;
            indexes_n[((n-1)*(N_X_-1) + m-1)*3+1] = n;
            indexes_l[((n-1)*(N_X_-1) + m-1)*3+1] = 2;
            values [((n-1)*(N_X_-1) + m-1)*3+2] = eval.lambda_3;
            indexes_m[((n-1)*(N_X_-1) + m-1)*3+2] = m;
            indexes_n[((n-1)*(N_X_-1) + m-1)*3+2] = n;
            indexes_l[((n-1)*(N_X_-1) + m-1)*3+2] = 3;
        }
    //print_evals(values, indexes_m, indexes_n, indexes_l);
    sort_evals(values, indexes_m, indexes_n, indexes_l);
    printf("\n\n");
    print_evals(values, indexes_m, indexes_n, indexes_l);

    return;
}

void arrays_init_for_second (double **u_out, double **u_in, double **u_curr, double **u_new,
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

void arrays_free_for_second (double **u_out, double **u_in, double **u_curr, double **u_new,
                  double **v_out, double **v_in, double **v_curr, double **v_new,
                  double **p_out, double **p_in, double **p_curr, double **p_new){

    free(*u_out); free(*u_in); free(*u_curr); free(*u_new);
    free(*v_out); free(*v_in); free(*v_curr); free(*v_new);
    free(*p_out); free(*p_in); free(*p_curr); free(*p_new);

    return;
}

void arrays_init_for_third (double **u_out, double **u_in, double **u_curr, double **u_new,
                            double **v_out, double **v_in, double **v_curr, double **v_new,
                            double **p_out, double **p_in, double **p_curr, double **p_new,
                            double **P_in, double **P_out, double **P,
                            double **C_in, double **C_out, double **C,
                            double **D_in, double **D_out, double **D,
                            int N_x, int N_y) {

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

    *P_in  = (double *) malloc (sizeof(double) * N_x * N_y);
    *P_out = (double *) malloc (sizeof(double) * N_x * N_y);
    *P     = (double *) malloc (sizeof(double) * N_x * N_y);

    *C_in  = (double *) malloc (sizeof(double) * N_x * N_y);
    *C_out = (double *) malloc (sizeof(double) * N_x * N_y);
    *C     = (double *) malloc (sizeof(double) * N_x * N_y);

    *D_in  = (double *) malloc (sizeof(double) * N_x * N_y);
    *D_out = (double *) malloc (sizeof(double) * N_x * N_y);
    *D     = (double *) malloc (sizeof(double) * N_x * N_y);
}
void arrays_free_for_third (double **u_out, double **u_in, double **u_curr, double **u_new,
                            double **v_out, double **v_in, double **v_curr, double **v_new,
                            double **p_out, double **p_in, double **p_curr, double **p_new,
                            double **P_in, double **P_out, double **P,
                            double **C_in, double **C_out, double **C,
                            double **D_in, double **D_out, double **D) {

    free(*u_out); free(*u_in); free(*u_curr); free(*u_new);
    free(*v_out); free(*v_in); free(*v_curr); free(*v_new);
    free(*p_out); free(*p_in); free(*p_curr); free(*p_new);

    free(*P_in); free(*P_out); free(*P);
    free(*C_in); free(*C_out); free(*C);
    free(*D_in); free(*D_out); free(*D);

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

void print_matr (double *a, int N, int M){
    FILE *fp;
    int i = 0, j = 0;

    fp = fopen("matrix","a");

    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j)
            fprintf(fp, "%+.5lf \t", a[i*M + j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n \n \n");
    fclose(fp);

    return;
}

void print_evals (double *eval, int *ind_m, int *ind_n, int *ind_l) {
    for (int n = 0; n < (N_Y_-1)*(N_X_-1)*3; ++n) {
            printf ("%d %d   %d \t %lf \n",
                        ind_m[n], ind_n[n], ind_l[n], eval[n]);

    }
    return;
}
void sort_evals (double *eval, int *ind_m, int *ind_n, int *ind_l) {
    //double min_value;
    double tmp;
    int tmp_;
    for(int i = 0 ; i < (N_Y_-1)*(N_X_-1)*3 - 1; i++) {
        for(int j = 0 ; j < (N_Y_-1)*(N_X_-1)*3 - i - 1 ; j++) {
            if(fabs(eval[j]) > fabs(eval[j+1])) {
                tmp = eval[j];
                eval[j] = eval[j+1] ;
                eval[j+1] = tmp;
                tmp_ = ind_m[j];
                ind_m[j] = ind_m[j+1];
                ind_m[j+1] = tmp_;
                tmp_ = ind_n[j];
                ind_n[j] = ind_n[j+1];
                ind_n[j+1] = tmp_;
                tmp_ = ind_l[j];
                ind_l[j] = ind_l[j+1];
                ind_l[j+1] = tmp_;
            }
        }
    }
}