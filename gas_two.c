//
// Created by vover on 3/5/18.
//

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "laspack/itersolv.h"
#include "gas_two.h"
#include "case.h"
#include "residuals.h"
#include "laspack/rtc.h"
#include "functions.h"
#include "gnuploting.h"


#define EPS 1e-8
#define EEPPSS 1e-15
#define MAX_ITER 2000


void param_dif (P_gas *p_d)
{
    p_d->Segm_X = 2*M_PI;
    p_d->Segm_Y = 2*M_PI;
    p_d->Segm_T = 10;
    p_d->mu = 0.1;
    p_d->p_ro = 10;
    p_d->omega = 1;

    return;
}

void param_she_step(P_she *p_s, P_gas p_d, int it_t, int it_sp)
{
    p_s->M_x = zero_spl_x * it_sp;
    p_s->M_y = zero_spl_y * it_sp;
    p_s->N = zero_spl_t * it_t;
    p_s->eta = 1;
    p_s->h_x = p_d.Segm_X / p_s->M_x;
    p_s->h_y = p_d.Segm_Y / p_s->M_y;
    p_s->tau = p_d.Segm_T / p_s->N;
    p_s->Dim = (p_s->M_x + 1) * (p_s->M_y + 1); //* (p_s->N + 1);

    return;
}

void param_t_const (T_const *t_c, P_she p_s, P_gas p_d)
{
    double thx_ = p_s.tau/p_s.h_x;
    double thy_ = p_s.tau/p_s.h_y;
    t_c->thx = thx_;
    t_c->thy = thy_;
    t_c->thx05 = thx_/2.;
    t_c->thy05 = thy_/2.;
    t_c->thx2 = 2.*thx_;
    t_c->thy2 = 2.*thy_;
    t_c->thx4 = 4.*thx_;
    t_c->thy4 = 4.*thy_;
    t_c->thx32 = 3.*thx_/2.;
    t_c->thy32 = 3.*thy_/2.;
    t_c->tau2 = 2. * p_s.tau;
    t_c->tau4 = 4. * p_s.tau;
    t_c->tau6 = 6. * p_s.tau;
    t_c->thxx8 = 8. * thx_/p_s.h_x;
    t_c->thxx6 = 6. * thx_/p_s.h_x;
    t_c->thyy8 = 8. * thy_/p_s.h_y;
    t_c->thyy6 = 6. * thy_/p_s.h_y;
    t_c->thxy = p_s.tau / (2. * p_s.h_x * p_s.h_y);
    t_c->Max = 3. * thx_ * p_d.p_ro;
    t_c->May = 3. * thy_ * p_d.p_ro;
    //t_c->Max = 0;
    //t_c->May = 0;
    return;
}

void param_MUM_const (MUM_const *MUM_c, P_she p_s, double GG, P_gas p_d)
{
    double tauhx = p_s.tau/(p_s.h_x * p_s.h_x);
    double tauhy = p_s.tau/(p_s.h_y * p_s.h_y);
    double M = p_d.mu * GG;
    MUM_c->MUM = p_d.mu * GG;
    MUM_c->MU8x = 8. * M * tauhx;
    MUM_c->MU6x = 6. * M * tauhx;
    MUM_c->MU8y = 8. * M * tauhy;
    MUM_c->MU6y = 6. * M * tauhy;
    MUM_c->MUv1 = 6. + M * (16. * tauhx + 12. * tauhy);
    MUM_c->MUv2 = 6. + M * (16. * tauhy + 12. * tauhx);
    return;
}

void param_MM_step (MM_step *MM_s, size_t mm, int n, int m)
{
    MM_s->mmg00 = mm;
    MM_s->mmv100 = mm + 1;
    MM_s->mmv200 = mm + 2;
    MM_s->mmgL0 = mm - 3;
    MM_s->mmv1L0 = mm - 2;
    MM_s->mmv2L0 = mm - 1;
    MM_s->mmgR0 = mm + 3;
    MM_s->mmv1R0 = mm + 4;
    MM_s->mmv2R0 = mm + 5;
    MM_s->mmg0L = (size_t) 3 * (m - n) + 1;
    MM_s->mmv10L = MM_s->mmg0L + 1;
    MM_s->mmv20L = MM_s->mmv10L + 1;
    MM_s->mmg0R = (size_t) 3 * (m + n) + 1;
    MM_s->mmv10R = MM_s->mmg0R + 1;
    MM_s->mmv20R = MM_s->mmv10R + 1;

    return;
}

int son (int i, int j, P_she *p_s) //  the status of the node
{
    if (i == 0 && j == 0) // x = 0; y = 0
        return 5;
    if (i == p_s->M_y && j == 0) // x = 0; y = b
        return 7;
    if (i == 0 && j == p_s->M_x) // x = a; y = 0;
        return 6;
    if (i == p_s->M_y && j == p_s->M_x) // x = a; y = b
        return 8;
    if (i > 0 && i < (p_s->M_y/2) && j == p_s->M_x/2) // Interpretation
        return 9;
    if (i == 0 && j > 0 && j < p_s->M_x/2) // x in; y = 0
        return 3;
    if (i == p_s->M_y && j > 0 && j < p_s->M_x) // x in; y = b
        return 4;
    if (i > 0 && i < p_s->M_y && j == 0) // x = 0; y in
        return 1;
    if (i > 0 && i < p_s->M_y && j == p_s->M_x) // x = a; y in
        return 2;
    if (i == 0 && j > p_s->M_x/2-1 && j < p_s->M_x)
        return 10;

    return 0; // internal mesh node
}

void print_vector(double* a, int n)
{
    int i;

    for (i = 0; i < n; ++i)
        printf ("%e \n", a[i]);
    printf("let's going on\n\n\n");

    return;
}

void Setka (int *st, P_she *p_s)
{
    int k = 0;

    for (int i = 0; i < p_s->M_y + 1; ++i) {
        for (int j = 0; j < p_s->M_x + 1; ++j) {
            st[k] = son(i, j, p_s);
 //           printf("%d \t", st[k]);
            k++;
        }
//        printf("\n");
    }
//    sleep(10);
    return;

}

void Sxema (double *G, double *V1, double *V2, int *st, P_she p_s, P_gas p_d)
{
    int k;

    QMatrix_L A;
    Vector b, x;

    T_const t_c;
    MUM_const m_c;
    MM_step m_s;

    double GG = 0;
    double tmp;
    size_t mm;
    int m;


    if (SMOOTH_SOLUTION == 1) first_fill(V1, V2, G, p_s, p_d.omega, u1, u2, g);
    else first_fill___ (V1, V2, G, p_s, p_d.omega);
/*
    print_vector(G, p_s.Dim);
    print_vector(V1, p_s.Dim);
    print_vector(V2, p_s.Dim);
*/

    param_t_const(&t_c, p_s, p_d);
    SetRTCAccuracy(EPS);


    for (k = 1; k < p_s.N + 1; ++k)
    {
        mm = 1;

        Q_Constr(&A, "A", (size_t) 3 * p_s.Dim, False, Rowws, Normal, True);
        V_Constr(&b, "b", (size_t) 3 * p_s.Dim, Normal, True);
        V_Constr(&x, "x", (size_t) 3 * p_s.Dim, Normal, True);

        GG = -1000000;
        for (m = 0; m < p_s.Dim; ++m)
        {
//            printf("%lf \n", GG);
            if (exp(-G[m]) > GG) GG = exp(-G[m]);
        }
        printf("%lf \n", GG);
        param_MUM_const(&m_c, p_s, GG, p_d);

        for (size_t i = 0; i < p_s.Dim; ++i)
        {
            V_SetCmp(&x, 3 * i + 1, G[i]);
            V_SetCmp(&x, 3 * i + 2, V1[i]);
            V_SetCmp(&x, 3 * i + 3, V2[i]);
        }


        for (m = 0; m < p_s.Dim; ++m)
        {
            param_MM_step(&m_s,mm,p_s.M_x+1, m);
            switch (st[m])
            {
                case 0:
                    mm = case0(&A, &b, t_c, m_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, p_d.p_ro, mm);
                    break;
                case 1:
                    mm = case1(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                case 2:
                    mm = case2(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                case 3:
                    mm = case3(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm, p_d.omega);
                    break;
                case 4:
                    mm = case4(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                case 5:
                    mm = case1(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                case 6:
                    mm = case2(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                case 7:
                    mm = case1(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                case 8:
                    mm = case2(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                case 9:
                    mm = case9(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                case 10:
                    mm = case10(&A, &b, t_c, m_s, k, V1, V2, G, m, p_s, p_d.mu, mm);
                    break;
                default:
                {
                    printf("FATAL ERROR: unknown type of node");
                    exit(1);
                }

            }
            ++mm;
        }

        CGSIter(&A, &x, &b, MAX_ITER, SSORPrecond, 1);

        for (size_t i = 0; i < p_s.Dim; ++i)
        {
            G[i] = V_GetCmp(&x, 3 * i + 1);
            V1[i] = V_GetCmp(&x, 3 * i + 1 + 1);
            V2[i] = V_GetCmp(&x, 3 * i + 2 + 1);
        }

        for (m = 0; m < p_s.Dim; ++m)
        {
            if (fabs(G[m])  < EEPPSS ) G[m] = 0;
            if (fabs(V1[m]) < EEPPSS ) V1[m] = 0;
            if (fabs(V2[m]) < EEPPSS ) V2[m] = 0;
        }
        sleep(1);


        Q_Destr(&A);
        V_Destr(&b);
        V_Destr(&x);

        if (SMOOTH_SOLUTION != 1)
            run_gnuplot(p_s, V1, V2, G, k);
    }



    return;
}

