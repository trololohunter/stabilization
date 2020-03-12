//
// Created by vover on 3/11/20.
//

#include <malloc.h>
#include "math.h"
#include "harmonic.h"
#include "gas_two.h"

void phi_mn_sin_cos (double *phi, int m, int n, P_she p_s){

    for (int i = 0; i < p_s.M_y + 1; ++i)
        for (int j = 0; j < p_s.M_x + 1; ++j)
            phi[i * (p_s.M_x + 1) + j] = sin (m * j * p_s.h_x / 2) * cos (n * i * p_s.h_y / 2);

    return;
}
void phi_mn_cos_sin (double *phi, int m, int n, P_she p_s){

    for (int i = 0; i < p_s.M_y + 1; ++i)
        for (int j = 0; j < p_s.M_x + 1; ++j)
            phi[i * (p_s.M_x + 1) + j] = cos (m * j * p_s.h_x / 2) * sin (n * i * p_s.h_y / 2);

    return;
}
void phi_mn_cos_cos (double *phi, int m, int n, P_she p_s){

    for (int i = 0; i < p_s.M_y + 1; ++i)
        for (int j = 0; j < p_s.M_x + 1; ++j)
            phi[i * (p_s.M_x + 1) + j] = cos (m * j * p_s.h_x / 2) * cos (n * i * p_s.h_y / 2);

    return;
}

double scalar_product (double *x, double *y, P_she p_s){

    double rez = 0;
    int i = 0;
    int j = 0;

    for (i = 0; i < p_s.M_x+1; ++i)
        rez += x[i] * y[i] * p_s.h_y * p_s.h_y / 2;

    for (j = 1; j < p_s.M_y; ++j)
        for (i = 0; i < p_s.M_x + 1; ++i)
            rez += x[j*(p_s.M_x + 1) + i] * y[j*(p_s.M_x + 1) + i] * p_s.h_y * p_s.h_x;

    j = p_s.M_y;
    for (i = 0; i < p_s.M_x + 1; ++i)
        rez += x[j*(p_s.M_x + 1) + i] * y[j*(p_s.M_x + 1) + i] * p_s.h_y * p_s.h_x / 4;

    return rez;
}
double coefficient_Cmn (double *u, P_she p_s, int m, int n){
    double *phi;
    double rez = 0;

    phi = (double *) malloc (p_s.Dim * sizeof(double));

    phi_mn_sin_cos(phi, m, n, p_s);

    rez = scalar_product(u, phi, p_s) / scalar_product(phi, phi, p_s);

    free(phi);

    return rez;
}
double coefficient_Dmn (double *v, P_she p_s, int m, int n){
    double *phi;
    double rez = 0;

    phi = (double *) malloc (p_s.Dim * sizeof(double));

    phi_mn_cos_sin(phi, m, n, p_s);

    rez = scalar_product(v, phi, p_s) / scalar_product(phi, phi, p_s);

    free(phi);

    return rez;
}
double coefficient_Pmn (double *p, P_she p_s, int m, int n){
    double *phi;
    double rez = 0;

    phi = (double *) malloc (p_s.Dim * sizeof(double));

    phi_mn_cos_cos(phi, m, n, p_s);

    rez = scalar_product(p, phi, p_s) / scalar_product(phi, phi, p_s);

    free(phi);

    return rez;
}
void eigenvalue_mn (eigenvalue *eval, P_gas p_g, int m, int n){

    double W = sqrt(- 9 * C_ro_0 * m * m * ro_0 * ro_0
                    - 9 * C_ro_0 * m * n * ro_0 * ro_0
                    + m*m*m*m*p_g.mu*p_g.mu
                    + 2*m*m*n*n*p_g.mu*p_g.mu
                    + n*n*n*n*p_g.mu*p_g.mu);
    eval->lambda_1 = - p_g.mu*(m*m + n*n)/2;
    eval->lambda_2 = W / (6 * ro_0) + 2 * eval->lambda_1 / 3;
    eval->lambda_3 = W / (6 * ro_0) - 2 * eval->lambda_1 / 3;

    return;
}
void eigenvector_mn (eigenvector *eval, P_gas p_g, int m, int n){

    double W = sqrt(- 9 * C_ro_0 * m * m * ro_0 * ro_0
                    - 9 * C_ro_0 * m * n * ro_0 * ro_0
                    + m*m*m*m*p_g.mu*p_g.mu
                    + 2*m*m*n*n*p_g.mu*p_g.mu
                    + n*n*n*n*p_g.mu*p_g.mu);
    eval->v1[0] = 0;
    eval->v1[1] = - n/m;
    eval->v1[2] = 1;
    eval->v2[0] = (m*m*p_g.mu + 2*n*n*p_g.mu + 2 * W) / (6 * C_ro_0 * m)
                 - n*p_g.mu*(-m*m*p_g.mu + m*n*p_g.mu - 2*n*n*p_g.mu - 2*W)
               / (6 * C_ro_0 * (2*m*m*p_g.mu - m*n*p_g.mu + n*n*p_g.mu + 2*W ));
    eval->v2[1] = (m*m*p_g.mu - m*n*p_g.mu + 2*n*n*p_g.mu + 2*W)
                  / (2*m*m*p_g.mu - n*m*p_g.mu + n*n*p_g.mu + 2*W);
    eval->v2[2] = 1;
    eval->v3[0] = (m*m*p_g.mu + 2*n*n*p_g.mu - 2 * W) / (6 * C_ro_0 * m)
                  - n*p_g.mu*(m*m*p_g.mu - m*n*p_g.mu + 2*n*n*p_g.mu - 2*W)
                    / (6 * C_ro_0 * (-2*m*m*p_g.mu + m*n*p_g.mu - n*n*p_g.mu + 2*W ));
    eval->v3[1] = (m*m*p_g.mu - m*n*p_g.mu + 2*n*n*p_g.mu - 2*W)
                  / (2*m*m*p_g.mu - n*m*p_g.mu + n*n*p_g.mu - 2*W);
    eval->v3[2] = 1;

    return;
}