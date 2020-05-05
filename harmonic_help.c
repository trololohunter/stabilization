//
// Created by vover on 3/11/20.
//

#include <malloc.h>
#include "math.h"
#include "harmonic_help.h"
#include "gas_two.h"
#include "functions.h"

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
void eigenvalue_mn (eigenvalue *eval, double mu, int m, int n){

    double W = sqrt(- 9 * C_ro_0 * m * m * ro_0 * ro_0
                    - 9 * C_ro_0 * m * n * ro_0 * ro_0
                    + m*m*m*m*mu*mu
                    + 2*m*m*n*n*mu*mu
                    + n*n*n*n*mu*mu);
    eval->lambda_1 = - mu*(m*m + n*n)/2;
    eval->lambda_2 = W / (6 * ro_0) + 2 * eval->lambda_1 / 3;
    eval->lambda_3 = W / (6 * ro_0) - 2 * eval->lambda_1 / 3;

    return;
}
void eigenvector_mn (eigenvector *evec, double mu, int m, int n){

    double W = sqrt(- 9 * C_ro_0 * m * m * ro_0 * ro_0
                    - 9 * C_ro_0 * m * n * ro_0 * ro_0
                    + m*m*m*m*mu*mu
                    + 2*m*m*n*n*mu*mu
                    + n*n*n*n*mu*mu);
    evec->v1[0] = 0;
    evec->v1[1] = -n/m;
    evec->v1[2] = 1;
    evec->v2[0] = (m*m*mu + 2*n*n*mu + 2 * W) / (6 * C_ro_0 * m)
                 - n*mu*(-m*m*mu + m*n*mu - 2*n*n*mu - 2*W)
               / (6 * C_ro_0 * (2*m*m*mu - m*n*mu + n*n*mu + 2*W ));
    evec->v2[1] = (m*m*mu - m*n*mu + 2*n*n*mu + 2*W)
                  / (2*m*m*mu - n*m*mu + n*n*mu + 2*W);
    evec->v2[2] = 1;
    evec->v3[0] = (m*m*mu + 2*n*n*mu - 2 * W) / (6 * C_ro_0 * m)
                  - n*mu*(m*m*mu - m*n*mu + 2*n*n*mu - 2*W)
                    / (6 * C_ro_0 * (-2*m*m*mu + m*n*mu - n*n*mu + 2*W ));
    evec->v3[1] = (m*m*mu - m*n*mu + 2*n*n*mu - 2*W)
                  / (2*m*m*mu - n*m*mu + n*n*mu - 2*W);
    evec->v3[2] = 1;

    return;
}

void fill_with_vector (double *G, double *V1, double *V2, P_she p_s, int m, int n, double vector[N_])
{
    int i, j;
    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j) {
            V1[i * (p_s.M_x + 1) + j] = u(vector[0], j * p_s.h_x, i * p_s.h_y, m, n);
            V2[i * (p_s.M_x + 1) + j] = v(vector[1], j * p_s.h_x, i * p_s.h_y, m, n);
             G[i * (p_s.M_x + 1) + j] = p(vector[2], j * p_s.h_x, i * p_s.h_y, m, n);
        }
    return;
}

void fill_with_vector_new (double *u_, double *v_, double *p_, int N_x, int N_y, double h_x, double h_y, int m, int n, double vector[N_]) {

    int i, j;

    for (i = 0; i < N_y; ++i)
        for (j = 0; j < N_x + 1; ++j) {
            u_[i * (N_x + 1) + j] = u(vector[0], j * h_x, (i + 1/2) * h_y, m, n);
        }

    for (i = 0; i < N_y + 1; ++i)
        for (j = 0; j < N_x; ++j) {
            v_[i * N_x + j] = v(vector[1], (j + 1/2) * h_x, i * h_y, m, n);
        }

    for (i = 0; i < N_y; ++i)
        for (j = 0; j < N_x; ++j) {
            p_[i * N_x + j] = p(vector[2], (j+1/2) * h_x, (i+1/2) * h_y, m, n);
        }


    return;
}

void phi_mn_sin_cos_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y) {

    int i, j;

    for (i = 0; i < N_y; ++i)
        for (j = 0; j < N_x + 1; ++j) {
            phi[i * (N_x + 1) + j] = sin(j * h_x * m/2)*cos((i + 1/2) * h_y*n/2);
        }

}

void phi_mn_cos_sin_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y){
    int i, j;

    for (i = 0; i < N_y + 1; ++i)
        for (j = 0; j < N_x; ++j) {
            phi[i * N_x + j] = cos((j + 1/2) * h_x * m/2)*sin(i * h_y*n/2);
        }

}
void phi_mn_cos_cos_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y){

    int i, j;
    for (i = 0; i < N_y; ++i)
        for (j = 0; j < N_x; ++j) {
            phi[i * N_x + j] = cos((j + 1/2) * h_x * m/2)*cos((i + 1/2) * h_y*n/2);
        }

}

double scalar_product_new (double *x, double *y, int size, double h_x, double h_y){
    double rez = 0;
    int i = 0;
    int j = 0;

    for (i = 0; i < size; ++i)
        rez += x[i] * y[i];

    return rez;
}

double coefficient_Cmn_new (double *u, int size, int N_x, int N_y, double h_x,  double h_y, int m, int n){
    double *phi;
    double rez = 0;

    phi = (double *) malloc (size * sizeof(double));

    phi_mn_sin_cos_new(phi, m, n, N_x, N_y, h_x, h_y);

    rez = scalar_product_new(u, phi, size, h_x, h_y) / scalar_product_new(phi, phi, size, h_x, h_y);

    free(phi);

    return rez;

}
double coefficient_Dmn_new (double *v, int size, int N_x, int N_y, double h_x, double h_y, int m, int n){
    double *phi;
    double rez = 0;

    phi = (double *) malloc (size * sizeof(double));

    phi_mn_cos_sin_new(phi, m, n, N_x, N_y, h_x, h_y);

    rez = scalar_product_new(v, phi, size, h_x, h_y) / scalar_product_new(phi, phi, size, h_x, h_y);

    free(phi);

    return rez;

}
double coefficient_Pmn_new (double *p, int size, int N_x, int N_y, double h_x, double h_y, int m, int n){
    double *phi;
    double rez = 0;

    phi = (double *) malloc (size * sizeof(double));

    phi_mn_cos_cos_new(phi, m, n, N_x, N_y, h_x, h_y);

    rez = scalar_product_new(p, phi, size, h_x, h_y) / scalar_product_new(phi, phi, size, h_x, h_y);

    free(phi);

    return rez;

}

double lambda_ (int i, double h) {
    return 2*sin(i*h/4)/h;
}

void fill_coef_for_harmonic_system_ij_in (coef_for_harmonic_system *CoefForHarmonicSystem,
                                          int i, int j, double h_x, double h_y, double rho_0,
                                          double C_rho, double tau, double mu){
    double lambda_ihx = lambda_(i,h_x);
    double lambda_jhy = lambda_(j,h_y);
    double coef43hyhx = 1 + (tau * mu / rho_0) * ((4/3) * lambda_jhy * lambda_jhy + lambda_ihx * lambda_ihx);
    double coef43hxhy = 1 + (tau * mu / rho_0) * ((4/3) * lambda_ihx * lambda_ihx + lambda_jhy * lambda_jhy);
    double coefhxhy = lambda_ihx * lambda_ihx + lambda_jhy * lambda_jhy;
    double W = coef43hyhx * coef43hxhy
               - (tau * tau) * (2/3) * C_rho * lambda_ihx * lambda_ihx * lambda_jhy * lambda_jhy
               + tau * tau * C_rho * lambda_jhy * lambda_jhy * coef43hxhy
                + tau * tau * C_rho * lambda_ihx * lambda_ihx * coef43hyhx
                - (1/9) * (tau * tau * mu * mu / (rho_0 * rho_0)) * lambda_ihx * lambda_ihx * lambda_jhy * lambda_jhy;

    CoefForHarmonicSystem->P[0] = (coef43hyhx * coef43hxhy
                                  - (1/9) * (tau * tau * mu * mu / (rho_0 * rho_0)) * lambda_ihx * lambda_ihx * lambda_jhy * lambda_jhy)/W;
    CoefForHarmonicSystem->P[1] = rho_0 * tau * lambda_ihx * coefhxhy / W;
    CoefForHarmonicSystem->P[2] = rho_0 * tau * lambda_jhy * coefhxhy / W;
    CoefForHarmonicSystem->C[0] = (tau / rho_0) * C_rho * lambda_ihx * coefhxhy / W;
    CoefForHarmonicSystem->C[1] = (coef43hyhx + tau * tau * C_rho * lambda_jhy * lambda_jhy) / W;
    CoefForHarmonicSystem->C[2] = (tau * tau * C_rho * lambda_ihx * coefhxhy
                                   + (1/3) * tau * mu / rho_0 * lambda_ihx * coefhxhy) / W;
    CoefForHarmonicSystem->D[0] = (tau / rho_0) * C_rho * lambda_jhy * coefhxhy / W;;
    CoefForHarmonicSystem->D[1] = (tau * tau * C_rho * lambda_ihx * coefhxhy
                                   + (1/3) * tau * mu / rho_0 * lambda_ihx * coefhxhy) / W;
    CoefForHarmonicSystem->D[2] = (coef43hxhy + tau * tau * C_rho * lambda_ihx * lambda_ihx) / W;

    return;
}

void coef_from_function (/*in*/ double *u, double *v, double *p,
                        /*out*/ double *C, double *D, double *P,
                                double h_x, double h_y, int N_x, int N_y) {
    int n = 0, m = 0;
    int u_size = (N_x + 1) * N_y;
    int v_size = (N_y + 1) * N_x;
    int p_size = N_x * N_y;

    C[0] = 0;
    D[0] = 0;
    P[0] = 0;

    for (n = 1; n < N_y; ++n)
        for (m = 0; m < N_x; ++m) {
            C[n*N_x+m] = coefficient_Cmn_new(u, u_size, N_x, N_y, h_x, h_y, m, n);
            D[n*N_x+m] = coefficient_Dmn_new(v, v_size, N_x, N_y, h_x, h_y, m, n);
            P[n*N_x+m] = coefficient_Pmn_new(p, p_size, N_x, N_y, h_x, h_y, m, n);
        }
    for (m = 1; m < N_x; ++m) {
        C[m] = coefficient_Cmn_new(u, u_size, N_x, N_y, h_x, h_y, m, 0);
        D[m] = coefficient_Dmn_new(v, v_size, N_x, N_y, h_x, h_y, m, 0);
        P[m] = coefficient_Pmn_new(p, p_size, N_x, N_y, h_x, h_y, m, 0);
    }

    return;
}
void funcion_from_coef (/*in*/ double *C, double *D, double *P,
                        /*out*/ double *u, double *v, double *p,
                               double h_x, double h_y, int N_x, int N_y){
    int i = 0, j = 0, m = 0, n = 0;

    // u from C
    for (i = 0; i < N_y; ++i) {
        u[i*(N_x+1)+0] = 0;
        u[i*(N_x+1)+N_x] = 0;
        for (j = 1; j < N_x; ++j) {
            u[i*(N_x+1) + j] = 0;
            for (n = 0; n < N_y; ++n) {
                for (m = 0; m < N_x; ++m) {
                    u[i*(N_x+1) + j] += C[n*N_x+m] * sin(j * h_x * m/2)*cos((i + 1/2) * h_y*n/2);
                }
            }
        }
    }

    //v from D
    for (j = 0; j < N_x; ++j) {
        v[j] = 0;
        v[N_y*N_x+j] = 0;
    }
    for (i = 1; i < N_y; ++i)
        for (j = 0; j < N_x; ++j) {
            v[i*N_x + j] = 0;
            for (n = 0; n < N_y; ++n) {
                for (m = 0; m < N_x; ++m) {
                    v[i*N_x + j] += D[n*N_x+m] * cos((j + 1/2) * h_x * m/2)*sin(i * h_y*n/2);
                }
            }
        }

    //p from P
    for (i = 0; i < N_y; ++i)
        for (j = 0; j < N_x; ++j) {
            p[i*N_x + j] = 0;
            for (n = 0; n < N_y; ++n) {
                for (m = 0; m < N_x; ++m) {
                    p[i*N_x + j] += P[n*N_x+m] * cos((j + 1/2) * h_x * m/2)*cos((i + 1/2) * h_y*n/2);
                }
            }
        }


    return;
}

void count_next_coef ( double *C_in, double *D_in, double *P_in,
                       double *C_out, double *D_out, double *P_out,
                       int N_x, int N_y, double h_x, double h_y,
                       double rho_0, double C_rho, double tau, double mu){
    int i = 0, j = 0;
    coef_for_harmonic_system CfHS;


    for (j = 0; j < N_y; ++j)
        for (i = 0; i < N_x; ++i){
            fill_coef_for_harmonic_system_ij_in(&CfHS, i, j, h_x, h_y, rho_0, C_rho, tau, mu);
            P_out[j*N_x +i] = CfHS.P[0] * P_in[j*N_x +i]
                              + CfHS.P[1] * C_in[j*N_x +i]
                              + CfHS.P[2] * D_in[j*N_x +i];
            C_out[j*N_x +i] = CfHS.C[0] * P_in[j*N_x +i]
                              + CfHS.C[1] * C_in[j*N_x +i]
                              + CfHS.C[2] * D_in[j*N_x +i];
            D_out[j*N_x +i] = CfHS.D[0] * P_in[j*N_x +i]
                              + CfHS.D[1] * C_in[j*N_x +i]
                              + CfHS.D[2] * D_in[j*N_x +i];
        }


    return;
}