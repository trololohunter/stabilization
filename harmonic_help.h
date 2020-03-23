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


//void constructor_eigenvector (eigenvector *a);
//void destructor_eigenvector (eigenvector *a);
void phi_mn_sin_cos (double *phi, int m, int n, P_she p_s);
void phi_mn_cos_sin (double *phi, int m, int n, P_she p_s);
void phi_mn_cos_cos (double *phi, int m, int n, P_she p_s);
double scalar_product (double *x, double *y, P_she p_s);
double coefficient_Cmn (double *u, P_she p_s, int m, int n);
double coefficient_Dmn (double *v, P_she p_s, int m, int n);
double coefficient_Pmn (double *p, P_she p_s, int m, int n);
void eigenvalue_mn (eigenvalue *eval, double mu, int m, int n);
void eigenvector_mn (eigenvector *evec, double mu, int m, int n);
void fill_with_vector (double *G, double *V1, double *V2, P_she p_s, int m, int n, double vector[N_]);
void fill_with_vector_new (double *u, double *v, double *p, int N_x, int N_y,double h_x, double h_y, int m, int n, double vector[N_]);

void phi_mn_sin_cos_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y);
void phi_mn_cos_sin_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y);
void phi_mn_cos_cos_new (double *phi, int m, int n, int N_x, int N_y, double h_x, double h_y);
double scalar_product_new (double *x, double *y, int size, double h_x, double h_y);
double coefficient_Cmn_new (double *u, int size, int N_x, int N_y, double h_x, double h_y, int m, int n);
double coefficient_Dmn_new (double *v, int size, int N_x, int N_y, double h_x, double h_y, int m, int n);
double coefficient_Pmn_new (double *p, int size, int N_x, int N_y, double h_x, double h_y, int m, int n);

#endif //UNTITLED2_HARMONIC_H
