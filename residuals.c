//
// Created by vover on 4/3/18.
//

#include <stdio.h>
#include <unistd.h>
#include "residuals.h"
#include "gas_two.h"

double residual_Ch_step (double *V1, double *V2, double *G, P_she p_s, int k,
                         func u1, func u2, func ro)
{
    double V1max = 0, V2max = 0, Gmax = 0;
    int i, j;

    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1max = (fabs(u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) > V1max) ? fabs(u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) : V1max;
            V2max = (fabs(u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) > V2max) ? fabs(u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) : V2max;
            Gmax = (fabs(ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]) > Gmax) ? fabs(ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]) : Gmax;
 /*           printf ("i = %d \t j = %d \n V1max = %e \t V2max = %e \t Gmax = %e \n", i, j,
                    fabs(u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]),
                    fabs(u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]),
                    fabs(ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]));
  */      }

//    printf ("V1max = %e \t V2max = %e \t Gmax = %e \n", V1max, V2max, Gmax);
    if (V1max > V2max && V1max > Gmax) return V1max;
    if (V2max > V1max && V2max > Gmax) return V2max;
    return Gmax;

}

double residual_Ch (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s,
                    func u1, func u2, func ro)
{
    double V1max = 0, V2max = 0, Gmax = 0;
    int i, j;

    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1max = (fabs(u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) > V1max) ? fabs(u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) : V1max;
            V2max = (fabs(u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) > V2max) ? fabs(u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) : V2max;
            Gmax = (fabs(ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]) > Gmax) ? fabs(ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]) : Gmax;
        }

    printf ("armadas V1max = %e \t V2max = %e \t Gmax = %e \n", V1max, V2max, Gmax);

    n_s->Gnorm = Gmax;
    n_s->V1norm = V1max;
    n_s->V2norm = V2max;

    if (V1max > V2max && V1max > Gmax) return V1max;
    if (V2max > V1max && V2max > Gmax) return V2max;
    return Gmax;
}

double residual_L2h_step (double *V1, double *V2, double *G, P_she p_s, int k,
                          func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);

}

double residual_L2h (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s,
                     func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += (ro(p_s.N * p_s.tau, j * p_s.h_x , i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    n_s->V1norm = sqrt(V1res);
    n_s->V2norm = sqrt(V2res);
    n_s->Gnorm  = sqrt(Gres);

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);
}



double residual_W12 (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += 2 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += 2 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += 2 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    n_s->V1norm = sqrt(V1res);
    n_s->V2norm = sqrt(V2res);
    n_s->Gnorm  = sqrt(Gres);

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);
}

double residual_Ch_Sokol (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s,
                    func u1, func u2, func ro)
{
    double V1max = 0, V2max = 0, Gmax = 0;
    int i, j;

    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1max = (fabs(u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) > V1max) ? fabs(u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) : V1max;
            V2max = (fabs(u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) > V2max) ? fabs(u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) : V2max;
            if ((j < p_s.M_x) && (i < p_s.M_y))
            Gmax = (fabs(ro(p_s.N * p_s.tau, j * p_s.h_x + p_s.h_x/2., i * p_s.h_y + p_s.h_y/2.) - G[i * (p_s.M_x + 1) + j]) > Gmax) ? fabs(ro(p_s.N * p_s.tau, j * p_s.h_x + p_s.h_x/2., i * p_s.h_y + p_s.h_y/2.) - G[i * (p_s.M_x + 1) + j]) : Gmax;
        }

    printf ("armadas V1max = %e \t V2max = %e \t Gmax = %e \n", V1max, V2max, Gmax);

    n_s->Gnorm = Gmax;
    n_s->V1norm = V1max;
    n_s->V2norm = V2max;

    if (V1max > V2max && V1max > Gmax) return V1max;
    if (V2max > V1max && V2max > Gmax) return V2max;
    return Gmax;
}

double residual_L2h_Sokol (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s,
                     func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += (ro(p_s.N * p_s.tau, (j-1) * p_s.h_x + p_s.h_x/2., (i-1) * p_s.h_y + p_s.h_y/2.) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, (j-1) * p_s.h_x + p_s.h_x/2., (i-1) * p_s.h_y + p_s.h_y/2.) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; i < p_s.M_x; ++i)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; i < p_s.M_x; ++i)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    n_s->V1norm = sqrt(V1res);
    n_s->V2norm = sqrt(V2res);
    n_s->Gnorm  = sqrt(Gres);

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);
}



double residual_W12_Sokol (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += 2 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += 2 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += 2 * (ro(p_s.N * p_s.tau, (j-1) * p_s.h_x + p_s.h_x/2., (i-1) * p_s.h_y + p_s.h_y/2.) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, (j-1) * p_s.h_x + p_s.h_x/2., (i-1) * p_s.h_y + p_s.h_y/2.) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    n_s->V1norm = sqrt(V1res);
    n_s->V2norm = sqrt(V2res);
    n_s->Gnorm  = sqrt(Gres);

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);
}