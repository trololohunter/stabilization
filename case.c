//
// Created by vover on 3/26/18.
//

#include "gas_two.h"
#include "laspack/itersolv.h"
#include "functions.h"
#include "case.h"

typedef struct
{
    double g00;
    double gL0;
    double gR0;
    double g0L;
    double g0R;
    double gLL;
    double gRL;
    double gLR;
    double gRR;
    double v100;
    double v1L0;
    double v1R0;
    double v10L;
    double v10R;
    double v1LL;
    double v1RL;
    double v1LR;
    double v1RR;
    double v200;
    double v2L0;
    double v2R0;
    double v20L;
    double v20R;
    double v2LL;
    double v2RL;
    double v2LR;
    double v2RR;

} gv1v2;

void param_gv1v2_0 (double *V1, double *V2, double *G, int m, int n, gv1v2 *v)
{
    v->g00 = G[m];
    v->gL0 = G[m-1];
    v->gR0 = G[m+1];
    v->g0L = G[m - n];
    v->g0R = G[m + n];
    v->gLL = G[m - n - 1];
    v->gRL = G[m - n + 1];
    v->gLR = G[m + n - 1];
    v->gRR = G[m + n + 1];
    v->v100 = V1[m];
    v->v1L0 = V1[m-1];
    v->v1R0 = V1[m+1];
    v->v10L = V1[m - n];
    v->v10R = V1[m + n];
    v->v1LL = V1[m - n - 1];
    v->v1RL = V1[m - n + 1];
    v->v1LR = V1[m + n - 1];
    v->v1RR = V1[m + n + 1];
    v->v200 = V2[m];
    v->v2L0 = V2[m-1];
    v->v2R0 = V2[m+1];
    v->v20L = V2[m - n];
    v->v20R = V2[m + n];
    v->v2LL = V2[m - n - 1];
    v->v2RL = V2[m - n + 1];
    v->v2LR = V2[m + n - 1];
    v->v2RR = V2[m + n + 1];

    return;
}

void param_gv1v2_1 (double *V1, double *V2, double *G, int m, int n, gv1v2 *v)
{
    v->g00 = G[m];
    v->gL0 = 0;
    v->gR0 = G[m+1];
    v->g0L = 0;
    v->g0R = 0;
    v->gLL = 0;
    v->gRL = 0;
    v->gLR = 0;
    v->gRR = 0;
    v->v100 = V1[m];
    v->v1L0 = 0;
    v->v1R0 = V1[m+1];
    v->v10L = 0;
    v->v10R = 0;
    v->v1LL = 0;
    v->v1RL = 0;
    v->v1LR = 0;
    v->v1RR = 0;
    v->v200 = V2[m];
    v->v2L0 = 0;
    v->v2R0 = V2[m+1];
    v->v20L = 0;
    v->v20R = 0;
    v->v2LL = 0;
    v->v2RL = 0;
    v->v2LR = 0;
    v->v2RR = 0;

    return;
}

void param_gv1v2_2 (double *V1, double *V2, double *G, int m, int n, gv1v2 *v)
{
    v->g00 = G[m];
    v->gL0 = G[m-1];
    v->gR0 = 0;
    v->g0L = 0;
    v->g0R = 0;
    v->gLL = 0;
    v->gRL = 0;
    v->gLR = 0;
    v->gRR = 0;
    v->v100 = V1[m];
    v->v1L0 = V1[m-1];
    v->v1R0 = 0;
    v->v10L = 0;
    v->v10R = 0;
    v->v1LL = 0;
    v->v1RL = 0;
    v->v1LR = 0;
    v->v1RR = 0;
    v->v200 = V2[m];
    v->v2L0 = V2[m-1];
    v->v2R0 = 0;
    v->v20L = 0;
    v->v20R = 0;
    v->v2LL = 0;
    v->v2RL = 0;
    v->v2LR = 0;
    v->v2RR = 0;

    return;
}

void param_gv1v2_3 (double *V1, double *V2, double *G, int m, int n, gv1v2 *v)
{
    v->g00 = G[m];
    v->gL0 = 0;
    v->gR0 = 0;
    v->g0L = 0;
    v->g0R = G[m + n];
    v->gLL = 0;
    v->gRL = 0;
    v->gLR = 0;
    v->gRR = 0;
    v->v100 = V1[m];
    v->v1L0 = 0;
    v->v1R0 = 0;
    v->v10L = 0;
    v->v10R = V1[m + n];
    v->v1LL = 0;
    v->v1RL = 0;
    v->v1LR = 0;
    v->v1RR = 0;
    v->v200 = V2[m];
    v->v2L0 = 0;
    v->v2R0 = 0;
    v->v20L = 0;
    v->v20R = V2[m + n];
    v->v2LL = 0;
    v->v2RL = 0;
    v->v2LR = 0;
    v->v2RR = 0;

    return;
}

void param_gv1v2_4 (double *V1, double *V2, double *G, int m, int n, gv1v2 *v)
{
    v->g00 = G[m];
    v->gL0 = 0;
    v->gR0 = 0;
    v->g0L = G[m - n];
    v->g0R = 0;
    v->gLL = 0;
    v->gRL = 0;
    v->gLR = 0;
    v->gRR = 0;
    v->v100 = V1[m];
    v->v1L0 = 0;
    v->v1R0 = 0;
    v->v10L = V1[m - n];
    v->v10R = 0;
    v->v1LL = 0;
    v->v1RL = 0;
    v->v1LR = 0;
    v->v1RR = 0;
    v->v200 = V2[m];
    v->v2L0 = 0;
    v->v2R0 = 0;
    v->v20L = V2[m - n];
    v->v20R = 0;
    v->v2LL = 0;
    v->v2RL = 0;
    v->v2LR = 0;
    v->v2RR = 0;

    return;
}


void first_fill(double *V1, double *V2, double *G, P_she p_s, double w, func u1, func u2, func ro) {
    int i, j;
    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1[i * (p_s.M_x + 1) + j] = u1(0, j * p_s.h_x, i * p_s.h_y);
            V2[i * (p_s.M_x + 1) + j] = u2(0, j * p_s.h_x, i * p_s.h_y);
            G[i * (p_s.M_x + 1) + j] = ro(0, j * p_s.h_x, i * p_s.h_y);
            if (i > 0 && i < (p_s.M_y/2) && j == p_s.M_x/2)
            {
                V1[i * (p_s.M_x + 1) + j] = 0;
                V2[i * (p_s.M_x + 1) + j] = 0;
            }
        }
    //for (j = 1; j < p_s.M_x/2; ++j)
    //V2[j] = w;
}

void first_fill___ (double *V1, double *V2, double *G, P_she p_s,
                    double w)
{
    int i, j;
    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1[i * (p_s.M_x + 1) + j] = 0;
            V2[i * (p_s.M_x + 1) + j] = 0;
            G[i * (p_s.M_x + 1) + j] = 0;
        }
    for (j = 1; j < p_s.M_x/2; ++j)
    V2[j] = w;

}

size_t case0 (QMatrix_L *A, Vector *B, T_const t_c, MUM_const m_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, double p_ro, size_t mm)
{
    double tmp;
    gv1v2 v;

    param_gv1v2_0(V1, V2, G, m, p_s.M_x + 1, &v);

    Q_SetLen(A, mm, 9);
    Q_SetEntry(A, mm, 0, mm, 4.);

    tmp = t_c.thx * (v.v100 + v.v1R0);
    Q_SetEntry(A, mm, 1, m_s.mmgR0, tmp);

    tmp = -t_c.thx * (v.v100 + v.v1L0);
    Q_SetEntry(A, mm, 2, m_s.mmgL0, tmp);

    tmp = t_c.thy * (v.v200 + v.v20R);
    Q_SetEntry(A, mm, 3, m_s.mmg0R, tmp);

    tmp = -t_c.thy * (v.v200 + v.v20L);
    Q_SetEntry(A, mm, 4, m_s.mmg0L, tmp);

    Q_SetEntry(A, mm, 5, m_s.mmv1R0, t_c.thx2);
    Q_SetEntry(A, mm, 6, m_s.mmv20R, t_c.thy2);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thx2);
    Q_SetEntry(A, mm, 8, m_s.mmv20L, -t_c.thy2);

    tmp = v.g00 * (4. + t_c.thx * (v.v1R0 - v.v1L0) + t_c.thy * (v.v20R - v.v20L))
            + t_c.tau4 * Func_0(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y);
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 7);
    Q_SetEntry(A, mm, 0, mm, m_c.MUv1);

    tmp = t_c.thx * (v.v100 + v.v1R0) - m_c.MU8x;
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, tmp);

    tmp = t_c.thy32 * (v.v200 + v.v20R) - m_c.MU6y;
    Q_SetEntry(A, mm, 2, m_s.mmv10R, tmp);

    tmp = t_c.Max;
    Q_SetEntry(A, mm, 3, m_s.mmgL0, -tmp);
    Q_SetEntry(A, mm, 4, m_s.mmgR0, tmp);

    tmp = -t_c.thx * (v.v100 + v.v1L0) - m_c.MU8x;
    Q_SetEntry(A, mm, 5, m_s.mmv1L0, tmp);

    tmp = -t_c.thy32 * (v.v200 + v.v20L) - m_c.MU6y;
    Q_SetEntry(A, mm, 6, m_s.mmv10L, tmp);

    tmp = v.v100 * (6. + t_c.thy32 *(v.v20R - v.v20L))
          - (m_c.MUM - mu * exp(-v.g00)) * (t_c.thxx8 * (v.v1R0 - 2 * v.v100 + v.v1L0) + t_c.thyy6 * (v.v10R - 2 * v.v100 + v.v10L))
          + t_c.thxy * (v.v2RR + v.v2LL - v.v2RL - v.v2LR) * mu * exp(-v.g00)
          + t_c.tau6 * Func_1(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y, p_ro, mu);
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 7);
    Q_SetEntry(A, mm, 0, mm, m_c.MUv2);

    tmp = t_c.thx32 * (v.v100 + v.v1R0) - m_c.MU6x;
    Q_SetEntry(A, mm, 1, m_s.mmv2R0, tmp);

    tmp = t_c.thy * (v.v200 + v.v20R) - m_c.MU8y;
    Q_SetEntry(A, mm, 2, m_s.mmv20R, tmp);

    tmp = t_c.May;
    Q_SetEntry(A, mm, 3, m_s.mmg0L, -tmp);
    Q_SetEntry(A, mm, 4, m_s.mmg0R, tmp);

    tmp = -t_c.thx32 * (v.v100 + v.v1L0) - m_c.MU6x;
    Q_SetEntry(A, mm, 5, m_s.mmv2L0, tmp);

    tmp = -t_c.thy * (v.v200 + v.v20L) - m_c.MU8y;
    Q_SetEntry(A, mm, 6, m_s.mmv20L, tmp);

    tmp = v.v200 * (6. + t_c.thx32 * (v.v1R0 - v.v1L0))
          - (m_c.MUM - mu * exp(-v.g00)) * (t_c.thxx6 * (v.v2R0 - 2 * v.v200 + v.v2L0) + t_c.thyy8 * (v.v20R - 2 * v.v200 + v.v20L))
          + t_c.thxy * (v.v1RR + v.v1LL - v.v1RL - v.v1LR) * mu * exp(-v.g00)
          + t_c.tau6 * Func_2(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y, p_ro, mu);
    V_SetCmp(B, mm, tmp);


    return mm;
}

size_t case1 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm)
{
    (void) mu;

    double tmp;
    gv1v2 v;

    param_gv1v2_1(V1, V2, G, m, p_s.M_x + 1, &v);

    Q_SetLen(A, mm, 3);
    Q_SetEntry(A, mm, 0, m_s.mmg00, 2.);
    tmp = t_c.thx * v.v1R0;
    Q_SetEntry(A, mm, 1, m_s.mmgR0, tmp);
    Q_SetEntry(A, mm, 2, m_s.mmv1R0, t_c.thx2);

    tmp = 2. * v.g00 + t_c.thx * v.g00 * v.v1R0
          + t_c.thx * (-2.5 * v.gR0 * v.v1R0 + 2. * G[m+2] * V1 [m+2] - 0.5 * G[m+3] * V1[m+3]
                         + (2 - v.g00) * (-2.5 * v.v1R0 + 2 * V1[m+2] - 0.5 * V1[m+3]))
          +  t_c.tau2 * Func_0(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv100, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv200, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);
    return mm;
}

size_t case2 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm)
{
    (void) mu;

    double tmp;
    int n = p_s.M_x + 1;
    gv1v2 v;

    param_gv1v2_2(V1, V2, G, m, n, &v);

    Q_SetLen(A, mm, 3);
    Q_SetEntry(A, mm, 0, m_s.mmg00, 2.);
    tmp = -t_c.thx * v.v1L0;
    Q_SetEntry(A, mm, 1, m_s.mmgL0, tmp);
    Q_SetEntry(A, mm, 2, m_s.mmv1L0, -t_c.thx2);

    tmp = 2. * v.g00 - t_c.thx * v.g00 * v.v1L0
          - t_c.thx * (-2.5 * v.gL0 * v.v1L0 + 2. * G[m-2] * V1 [m-2] - 0.5 * G[m-3] * V1[m-3]
                         + (2 - v.g00) * (-2.5 * v.v1L0 + 2. * V1[m-2] - 0.5 * V1[m-3]))
          +  t_c.tau2 * Func_0(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv100, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv200, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);

    return mm;
}


size_t case3 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm, double w)
{
    (void) mu;

    double tmp;
    int n = p_s.M_x + 1;
    gv1v2 v;

    param_gv1v2_3(V1, V2, G, m, n, &v);

    if (SMOOTH_SOLUTION == 1) {
        Q_SetLen(A, mm, 3);
        Q_SetEntry(A, mm, 0, mm, 2.);
        tmp = t_c.thy * v.v20R;
        Q_SetEntry(A, mm, 1, m_s.mmg0R, tmp);
        Q_SetEntry(A, mm, 2, m_s.mmv20R, t_c.thy2);
    }
    else {
        Q_SetLen(A, mm, 4);
        tmp = 2. - t_c.thy * v.v200;
        Q_SetEntry(A, mm, 0, m_s.mmg00, tmp);
        tmp = t_c.thy * v.v20R;
        Q_SetEntry(A, mm, 1, m_s.mmg0R, tmp);
        Q_SetEntry(A, mm, 2, m_s.mmv20R, t_c.thy2);
        Q_SetEntry(A, mm, 3, m_s.mmv200, -t_c.thy2);
    }

    tmp = 2. * v.g00 + t_c.thy * v.g00 * (v.v20R - v.v200)
          + t_c.thy * (v.g00 * v.v200 -2.5 * v.g0R * v.v20R + 2. * G[m+2*n] * V2 [m+2*n] - 0.5 * G[m+3*n] * V2[m+3*n]
                       + (2 - v.g00) * (v.v200 -2.5 * v.v20R + 2 * V2[m+2*n] - 0.5 * V2[m+3*n]))
          +  t_c.tau2 * Func_0(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv100, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv200, 1.);
    if (SMOOTH_SOLUTION == 1) {
        tmp = 0; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
        (void) w;
    }
    else tmp = w;
    V_SetCmp(B, mm, tmp);

    return mm;
}

size_t case4 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm)
{
    (void) mu;

    double tmp;
    int n = p_s.M_x + 1;
    gv1v2 v;

    param_gv1v2_4(V1, V2, G, m, n, &v);

    Q_SetLen(A, mm, 3);
    Q_SetEntry(A, mm, 0, m_s.mmg00, 2.);
    tmp = -t_c.thy * v.v20L;
    Q_SetEntry(A, mm, 1, m_s.mmg0L, tmp);
    Q_SetEntry(A, mm, 2, m_s.mmv20L, -t_c.thy2);

    tmp = 2. * v.g00 - t_c.thy * v.g00 * v.v20L
          - t_c.thy * (-2.5 * v.g0L * v.v20L + 2. * G[m-2*n] * V2 [m-2*n] - 0.5 * G[m-3*n] * V2[m-3*n]
                       + (2 - v.g00) * (-2.5 * v.v20L + 2 * V2[m-2*n] - 0.5 * V2[m-3*n]))
          +  t_c.tau2 * Func_0(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv100, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv200, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);

    return mm;
}

size_t case9 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm)
{
    (void) mu;

    double tmp;
    int n = p_s.M_x + 1;
    gv1v2 v;

    param_gv1v2_0(V1, V2, G, m, n, &v);

    Q_SetLen(A, mm, 9);
    Q_SetEntry(A, mm, 0, mm, 4.);

    tmp = t_c.thx * (v.v100 + v.v1R0);
    Q_SetEntry(A, mm, 1, m_s.mmgR0, tmp);

    tmp = -t_c.thx * (v.v100 + v.v1L0);
    Q_SetEntry(A, mm, 2, m_s.mmgL0, tmp);

    tmp = t_c.thy * (v.v200 + v.v20R);
    Q_SetEntry(A, mm, 3, m_s.mmg0R, tmp);

    tmp = -t_c.thy * (v.v200 + v.v20L);
    Q_SetEntry(A, mm, 4, m_s.mmg0L, tmp);

    Q_SetEntry(A, mm, 5, m_s.mmv1R0, t_c.thx2);
    Q_SetEntry(A, mm, 6, m_s.mmv20R, t_c.thy2);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thx2);
    Q_SetEntry(A, mm, 8, m_s.mmv20L, -t_c.thy2);

    tmp = v.g00 * (4. + t_c.thx * (v.v1R0 - v.v1L0) + t_c.thy * (v.v20R - v.v20L))
            + t_c.tau4 * Func_0(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y);
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv100, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv200, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);


    return mm;
}

size_t case10 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm)
{
    (void) mu;

    double tmp;
    int n = p_s.M_x + 1;
    gv1v2 v;

    param_gv1v2_3(V1, V2, G, m, n, &v);

    if (SMOOTH_SOLUTION == 1) {
        Q_SetLen(A, mm, 3);
        Q_SetEntry(A, mm, 0, mm, 2.);
        tmp = t_c.thy * v.v20R;
        Q_SetEntry(A, mm, 1, m_s.mmg0R, tmp);
        Q_SetEntry(A, mm, 2, m_s.mmv20R, t_c.thy2);
    }
    else {
        Q_SetLen(A, mm, 4);
        tmp = 2. - t_c.thy * v.v200;
        Q_SetEntry(A, mm, 0, m_s.mmg00, tmp);
        tmp = t_c.thy * v.v20R;
        Q_SetEntry(A, mm, 1, m_s.mmg0R, tmp);
        Q_SetEntry(A, mm, 2, m_s.mmv20R, t_c.thy2);
        Q_SetEntry(A, mm, 3, m_s.mmv200, -t_c.thy2);
    }

    tmp = 2. * v.g00 + t_c.thy * v.g00 * (v.v20R - v.v200)
          + t_c.thy * (v.g00 * v.v200 - 2.5 * v.g0R * v.v20R
                       + 2. * G[m+2*n] * V2 [m+2*n] - 0.5 * G[m+3*n] * V2[m+3*n]
                       + (2 - v.g00) * (v.v200 -2.5 * v.v20R
                                        + 2 * V2[m+2*n] - 0.5 * V2[m+3*n]))
          +  t_c.tau2 * Func_0(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y);
    V_SetCmp(B, mm, tmp);

    mm++;
    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, m_s.mmv100, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(B, mm, tmp);

    mm++;
    if (SMOOTH_SOLUTION == 1) {
        Q_SetLen(A, mm, 1);
        Q_SetEntry(A, mm, 0, m_s.mmv200, 1.);
        tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
        V_SetCmp(B, mm, tmp);
    }
    else {
        Q_SetLen(A, mm, 2);
        Q_SetEntry(A, mm, 0, m_s.mmv200, -1.);
        Q_SetEntry(A, mm, 1, m_s.mmv20R, 1.);
        tmp = 0; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
        V_SetCmp(B, mm, tmp);
    }

    return mm;
}