//
// Created by vover on 3/11/20.
//

#include "gas_two.h"

#ifndef UNTITLED2_HARMONIC_H
#define UNTITLED2_HARMONIC_H

#define N_ 3
#define ro_0 0.01
#define C_ro_0 1

typedef struct {
    double lambda_1;
    double lambda_2;
    double lambda_3;
} eigenvalue;

typedef struct {
    double v1[N_];
    double v2[N_];
    double v3[N_];
} eigenvector;

typedef struct {
    double P[N_];
    double C[N_];
    double D[N_];
} coef_for_harmonic_system;

//void constructor_eigenvector (eigenvector *a);
//void destructor_eigenvector (eigenvector *a);

void eigenvalue_mn (eigenvalue *eval, double mu, int m, int n);
void eigenvector_mn (eigenvector *evec, double mu, int m, int n);
void eigenvector_mn_t (eigenvector *evec, double mu, int m, int n);

void phi_mn_sin_cos (double *phi, int m, int n, P_she p_s);
void phi_mn_cos_sin (double *phi, int m, int n, P_she p_s);
void phi_mn_cos_cos (double *phi, int m, int n, P_she p_s);
double scalar_product (double *x, double *y, P_she p_s);
double coefficient_Cmn (double *u, P_she p_s, int m, int n);
double coefficient_Dmn (double *v, P_she p_s, int m, int n);
double coefficient_Pmn (double *p, P_she p_s, int m, int n);
void fill_with_vector (double *G, double *V1, double *V2, P_she p_s, int m, int n, double vector[N_]);

void fill_with_vector_new (double *u, double *v, double *p, int N_x, int N_y,double h_x, double h_y, int m, int n, double vector[N_]);
void phi_mn_sin_cos_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y);
void phi_mn_cos_sin_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y);
void phi_mn_cos_cos_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y);
double scalar_product_new (double *x, double *y, int size, double h_x, double h_y);
double coefficient_Cmn_new (double *u, int size, int N_x, int N_y, double h_x, double h_y, int m, int n);
double coefficient_Dmn_new (double *v, int size, int N_x, int N_y, double h_x, double h_y, int m, int n);
double coefficient_Pmn_new (double *p, int size, int N_x, int N_y, double h_x, double h_y, int m, int n);

double lambda_ (int i, double h);
void fill_coef_for_harmonic_system_ij_in (coef_for_harmonic_system *CoefForHarmonicSystem, int i, int j, double h_x, double h_y, double rho_0, double C_rho, double tau, double mu);

void coef_from_function (/*in*/ double *u, double *v, double *p,
                        /*out*/ double *C, double *D, double *P,
                                double h_x, double h_y, int N_x, int N_y);
void funcion_from_coef (/*in*/ double *C, double *D, double *P,
                        /*out*/ double *u, double *v, double *p,
                               double h_x, double h_y, int N_x, int N_y, double mu, double t);

void count_next_coef ( double *C_in, double *D_in, double *P_in,
                    double *C_out, double *D_out, double *P_out,
                    int N_x, int N_y, double h_x, double h_y,
                    double rho_0, double C_rho, double tau, double mu);

#endif //UNTITLED2_HARMONIC_H
