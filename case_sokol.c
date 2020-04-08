//
// Created by vover on 6/4/18.
//

#include "gas_two.h"
#include "laspack/itersolv.h"
#include "functions.h"
#include "case_sokol.h"

void first_fill_sokol (double *V1, double *V2, double *G, P_she p_s,
                       double w, func u1, func u2, func ro)
{
    int i, j;
    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1[i * (p_s.M_x + 1) + j] = u1(0, j * p_s.h_x, i * p_s.h_y);
            V2[i * (p_s.M_x + 1) + j] = u2(0, j * p_s.h_x, i * p_s.h_y);
            if ((j < p_s.M_x) && (i < p_s.M_y))
            G[i * (p_s.M_x + 1) + j] = ro(0, j * p_s.h_x + p_s.h_x/2., i * p_s.h_y + p_s.h_y/2.);
            if (i > 0 && i < (p_s.M_y/2) && j == p_s.M_x/2)
            {
                V1[i * (p_s.M_x + 1) + j] = 0;
                V2[i * (p_s.M_x + 1) + j] = 0;
            }
        }
}
void first_fill_sokol___ (double *V1, double *V2, double *G, P_she p_s,
                          double w)
{
    int i, j;
    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1[i * (p_s.M_x + 1) + j] = 0;
            V2[i * (p_s.M_x + 1) + j] = 0;
            if ((j < p_s.M_x) && (i < p_s.M_y))
            G[i * (p_s.M_x + 1) + j] = 0;
        }
    for (j = 1; j < p_s.M_x/2; ++j)
        V2[j] = w;
}
void case_sokol_H (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                   double *V1, double *V2, double *G,
                   int m, P_she p_s, double mu, size_t mm)
{
    size_t mmm = mm / (p_s.M_x) * (p_s.M_x + 1) + mm % (p_s.M_x);
    int n = p_s.M_x + 1;

    double V1s2 = (V1[mmm + n] + V1[mmm]) / 2.;
    double V2s1 = (V2[mmm + 1] + V2[mmm]) / 2.;
    double V1s2R = (V1[mmm + 1 + n] + V1[mmm + 1]) / 2.;
    double V2s1T = (V2[mmm + 1 + n] + V2[mmm + n]) / 2.;

    double tmp = 1;

    Q_SetLen(&H, mm, 3);
    if (V1s2 < 0)
    {
        Q_SetEntry(&H, mm, 1, mm+1, t_c.thx * V1s2R);
        tmp += -t_c.thx * V1s2;
    }
    else
    {
        Q_SetEntry(&H, mm, 1, mm-1, -t_c.thx * V1s2);
        tmp += t_c.thx * V1s2R;
    }

    if (V2s1 < 0)
    {
        Q_SetEntry(&H, mm, 2, mm+p_s.M_x, t_c.thy * V2s1T);
        tmp += -t_c.thy * V2s1;
    }
    else
    {
        Q_SetEntry(&H, mm, 2, mm-p_s.M_x, -t_c.thy * V2s1);
        tmp += t_c.thy * V2s1T;
    }
    Q_SetEntry(&H, mm, 0, mm, tmp);

    tmp = G[mm] + Func_0(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x + p_s.h_x/2., ((int)(m / (p_s.M_x + 1))) * p_s.h_y + p_s.h_y/2.);
    V_SetCmp(&b_H, mm, tmp);
}
void case0_sokol_V (QMatrix_L *U, Vector *b_U, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm, double p_ro)
{
    size_t mmm = m / (p_s.M_x+1) * (p_s.M_x) + m % (p_s.M_x+1);
    double Htilda = (Hnew[mmm] + Hnew[mmm-1] + Hnew[mmm-p_s.M_x] + Hnew[mmm-p_s.M_x-1])/4.;
    double tmp;
    double tauhx = p_s.tau/(p_s.h_x * p_s.h_x);
    double tauhy = p_s.tau/(p_s.h_y * p_s.h_y);
    double tauhxy = p_s.tau / (2. * p_s.h_x * p_s.h_y);
    int n = p_s.M_x + 1;
    double gamma = GAMMA;

    Q_SetLen(&U, mm, 7);

    if (V1[m] < 0) {

        tmp = 2 * tauhy * mu + 8. / 3. * tauhx * mu
              + Htilda - (Hnew[mmm] + Hnew[mmm - p_s.M_x]) * V1[m] * t_c.thx /4.
                - (Hnew[mmm-1] + Hnew[mmm-1-p_s.M_x]) * V1[m-1] * t_c.thx /4.;
        Q_SetEntry(&U, mm, 0, mm, tmp);

        tmp = (Hnew[mmm] + Hnew[mmm-p_s.M_x]) * V1[m] * t_c.thx /4.
                + (Hnew[mmm+1] + Hnew[mmm+1-p_s.M_x]) * V1[m+1] * t_c.thx /4.
                - tauhx * mu * 4./3.;
        Q_SetEntry(&U, mm, 1, mm+2, tmp);

        tmp = -tauhx * mu * 4./3.;
        Q_SetEntry(&U, mm, 2, mm-2, tmp);

        tmp = -tauhy * mu;
        Q_SetEntry(&U, mm, 3, 2*(mm-n) + 1, tmp);
        Q_SetEntry(&U, mm, 4, 2*(mm+n) + 1, tmp);

        tmp = -(Hnew[mmm+p_s.M_x] + Hnew[mmm-1+p_s.M_x]) * V1[m+n] * tauhy/4.
                - (Hnew[mmm] + Hnew[mmm-1]) * V1[m] * tauhy/4.;
        Q_SetEntry(&U, mm, 5, mm + 1, tmp);

        tmp = +(Hnew[mmm] + Hnew[mmm-1]) * V1[m] * tauhy/4.
              + (Hnew[mmm-p_s.M_x] + Hnew[mmm-1-p_s.M_x]) * V1[m-n] * tauhy/4.;
        Q_SetEntry(&U, mm, 6, 2*(mm+n)+2, tmp);

        tmp = tauhxy / 6. * mu * (V2[m+1+n] + V2[m-m-1] - V2[m-n+1] - V2[m+n-1])
                + p_s.tau * Htilda * Func_1(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y, p_ro, mu)
                - gamma/(gamma - 1) * Htilda *
                  (pow(Hnew[mmm] + Hnew[mmm-1],gamma-1)
                   - pow(Hnew[mmm-p_s.M_x] + Hnew[mmm-p_s.M_x-1],gamma-1))
                  * 1/p_s.h_x * pow(2, 1-gamma)
                + (G[mmm] + G[mmm-1] + G[mmm - p_s.M_x] + G[mmm-p_s.M_x-1])/4. * V1[m];
        V_SetCmp(b_U, mm, tmp);
    }
    else {
        tmp = 2 * tauhy * mu + 8. / 3. * tauhx * mu
              + Htilda + (Hnew[mmm] + Hnew[mmm - p_s.M_x]) * V1[m+1] * t_c.thx /4.
              + (Hnew[mmm-1] + Hnew[mmm-1-p_s.M_x]) * V1[m] * t_c.thx /4.;
        Q_SetEntry(&U, mm, 0, mm, tmp);

        tmp = -(Hnew[mmm-1] + Hnew[mmm-p_s.M_x-1]) * V1[m] * t_c.thx /4.
              - (Hnew[mmm-2] + Hnew[mmm-2-p_s.M_x]) * V1[m-1] * t_c.thx /4.
              - tauhx * mu * 4./3.;
        Q_SetEntry(&U, mm, 1, mm-2, tmp);

        tmp = -tauhx * mu * 4./3.;
        Q_SetEntry(&U, mm, 2, mm+2, tmp);

        tmp = -tauhy * mu;
        Q_SetEntry(&U, mm, 3, 2*(mm-n) + 1, tmp);
        Q_SetEntry(&U, mm, 4, 2*(mm+n) + 1, tmp);

        tmp = +(Hnew[mmm] + Hnew[mmm-1]) * V1[mm+n] * tauhy/4.
              + (Hnew[mmm-p_s.M_x] + Hnew[mmm-1-p_s.M_x]) * V1[m] * tauhy/4.;
        Q_SetEntry(&U, mm, 5, mm + 1, tmp);

        tmp = -(Hnew[mmm-1-p_s.M_x] + Hnew[mmm-1]) * V1[m] * tauhy/4.
              - (Hnew[mmm-2*p_s.M_x] + Hnew[mmm-1-2*p_s.M_x]) * V1[m-n] * tauhy/4.;
        Q_SetEntry(&U, mm, 6, 2*(mm-n)+2, tmp);

        tmp = tauhxy / 6. * mu * (V2[m+1+n] + V2[m-m-1] - V2[m-n+1] - V2[m+n-1])
              + p_s.tau * Htilda * Func_1(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y, p_ro, mu)
              - gamma/(gamma - 1) * Htilda *
                (pow(Hnew[mmm] + Hnew[mmm-1],gamma-1)
                 - pow(Hnew[mmm-p_s.M_x] + Hnew[mmm-p_s.M_x-1],gamma-1))
                * 1/p_s.h_x * pow(2, 1-gamma)
              + (G[mmm] + G[mmm-1] + G[mmm - p_s.M_x] + G[mmm-p_s.M_x-1])/4. * V1[m];
        V_SetCmp(b_U, mm, tmp);
    }
    ++mm;

    Q_SetLen(&U, mm, 7);
    if (V2[mm] < 0)
    {
        tmp = 2 * tauhx * mu + 8. / 3. * tauhy * mu
              + Htilda - (Hnew[mmm] + Hnew[mmm-1]) * V2[m] * tauhy / 4.
                - (Hnew[mmm-p_s.M_x] + Hnew[mmm-p_s.M_x-1]) * V2[m-n] * tauhy / 4.;
        Q_SetEntry(&U, mm, 0, mm, tmp);

        tmp = (Hnew[mmm+p_s.M_x] + Hnew[mmm+p_s.M_x-1]) * V2[m] * tauhy / 4.
                + (Hnew[mmm] + Hnew[mmm-1]) * V2[m+n] * tauhy / 4.
                - tauhy * mu * 4./3.;
        Q_SetEntry(&U, mm, 1, 2*(mm-1+n)+2, tmp);

        tmp = -tauhy * mu * 4./3.;
        Q_SetEntry(&U, mm, 2, 2*(mm-1-n)+2, tmp);

        tmp = -tauhx * mu;
        Q_SetEntry(&U, mm, 3, mm+2, tmp);
        Q_SetEntry(&U, mm, 4, mm-2, tmp);

        tmp = -(Hnew[mmm] + Hnew[mmm-p_s.M_x]) * V2[m] * tauhx / 4.
                - (Hnew[mmm-1] + Hnew[mmm-1-p_s.M_x]) * V2[m-1] * tauhx / 4.;
        Q_SetEntry(&U, mm, 5, mm - 1, tmp);

        tmp = (Hnew[mmm+1] + Hnew[mmm-p_s.M_x+1]) * V2[m+1] * tauhx / 4.
              + (Hnew[mmm] + Hnew[mmm-p_s.M_x]) * V2[m] * tauhx / 4.;
        Q_SetEntry(&U, mm, 6, mm + 1, tmp);

        tmp = tauhxy / 6. * mu * (V1[m+1+n] + V1[m-m-1] - V1[m-n+1] - V1[m+n-1])
              + p_s.tau * Htilda * Func_2(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y, p_ro, mu)
              - gamma/(gamma - 1) * Htilda *
                (pow(Hnew[mmm] + Hnew[mmm-1],gamma-1)
                 - pow(Hnew[mmm-p_s.M_x] + Hnew[mmm-p_s.M_x-1],gamma-1))
                * 1/p_s.h_x * pow(2, 1-gamma)
              + (G[mmm] + G[mmm-1] + G[mmm - p_s.M_x] + G[mmm-p_s.M_x-1])/4. * V2[m];
        V_SetCmp(b_U, mm, tmp);

    }
    else
    {
        tmp = 2 * tauhx * mu + 8. / 3. * tauhy * mu
              + Htilda + (Hnew[mmm] + Hnew[mmm-1]) * V2[m+n] * tauhy / 4.
              + (Hnew[mmm-p_s.M_x] + Hnew[mmm-p_s.M_x-1]) * V2[m] * tauhy / 4.;
        Q_SetEntry(&U, mm, 0, mm, tmp);

        tmp = -(Hnew[mmm-p_s.M_x] + Hnew[mmm-p_s.M_x-1]) * V2[m] * tauhy / 4.
              - (Hnew[mmm-2*p_s.M_x] + Hnew[mmm-2*p_s.M_x-1]) * V2[m-n] * tauhy / 4.
              - tauhy * mu * 4./3.;
        Q_SetEntry(&U, mm, 1, 2*(mm-1-n)+2, tmp);

        tmp = -tauhy * mu * 4./3.;
        Q_SetEntry(&U, mm, 2, 2*(mm-1+n)+2, tmp);

        tmp = -tauhx * mu;
        Q_SetEntry(&U, mm, 3, mm+2, tmp);
        Q_SetEntry(&U, mm, 4, mm-2, tmp);

        tmp = -(Hnew[mmm-1] + Hnew[mmm-p_s.M_x-1]) * V2[m] * tauhx / 4.
              - (Hnew[mmm-1-2*p_s.M_x] + Hnew[mmm-1-2*p_s.M_x]) * V2[m-1] * tauhx / 4.;
        Q_SetEntry(&U, mm, 5, mm - 3, tmp);

        tmp = (Hnew[mmm] + Hnew[mmm-p_s.M_x]) * V2[m+1] * tauhx / 4.
              + (Hnew[mmm-1] + Hnew[mmm-p_s.M_x-1]) * V2[m] * tauhx / 4.;
        Q_SetEntry(&U, mm, 6, mm - 1, tmp);

        tmp = tauhxy / 6. * mu * (V1[m+1+n] + V1[m-m-1] - V1[m-n+1] - V1[m+n-1])
              + p_s.tau * Htilda * Func_2(k * p_s.tau, m % (p_s.M_x + 1) * p_s.h_x, ((int)(m / (p_s.M_x + 1))) * p_s.h_y, p_ro, mu)
              - gamma/(gamma - 1) * Htilda *
                (pow(Hnew[mmm] + Hnew[mmm-1],gamma-1)
                 - pow(Hnew[mmm-p_s.M_x] + Hnew[mmm-p_s.M_x-1],gamma-1))
                * 1/p_s.h_x * pow(2, 1-gamma)
              + (G[mmm] + G[mmm-1] + G[mmm - p_s.M_x] + G[mmm-p_s.M_x-1])/4. * V2[m];
        V_SetCmp(b_U, mm, tmp);
    }
}
void case1_sokol_V (QMatrix_L *U, Vector *b_U, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm)
{
    double tmp;

    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);

    mm++;
    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);
}
void case2_sokol_V (QMatrix_L *U, Vector *b_U, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm)
{
    double tmp;

    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);

    mm++;
    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);
}
void case3_sokol_V (QMatrix_L *U, Vector *b_U, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm, double w)
{
    double tmp;

    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);

    mm++;
    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    if (SMOOTH_SOLUTION == 1) {
        tmp = 0; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
        (void) w;
    }
    else tmp = w;
    V_SetCmp(b_U, mm, tmp);
}
void case4_sokol_V (QMatrix_L *U, Vector *b_U, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm)
{
    double tmp;

    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);

    mm++;
    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);
}
void case9_sokol_V (QMatrix_L *U, Vector *b_U, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm)
{
    double tmp;

    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);

    mm++;
    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);
}
void case10_sokol_V (QMatrix_L *U, Vector *b_U, T_const t_c, int k,
                     double *V1, double *V2, double *G, double *Hnew,
                     int m, P_she p_s, double mu, size_t mm)
{
    double tmp;

    Q_SetLen(U, mm, 1);
    Q_SetEntry(U, mm, 0, mm, 1.);
    tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
    V_SetCmp(b_U, mm, tmp);

    mm++;
    if (SMOOTH_SOLUTION == 1) {
        Q_SetLen(U, mm, 1);
        Q_SetEntry(U, mm, 0, mm, 1.);
        tmp = 0.; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
        V_SetCmp(b_U, mm, tmp);
    }
    else {
        Q_SetLen(U, mm, 2);
        Q_SetEntry(U, mm, 0, mm, -1.);
        Q_SetEntry(U, mm, 1, mm + p_s.M_x + 1, 1.);
        tmp = 0; // +  t_c.tau6 * FUNC_2(tt, xx, yy);
        V_SetCmp(b_U, mm, tmp);
    }
}
