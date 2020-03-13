//
// Created by vover on 3/11/20.
//

#include "gas_two.h"

#ifndef UNTITLED2_HARMONIC_H
#define UNTITLED2_HARMONIC_H

#define N_ 3
#define ro_0 0.1
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
void eigenvalue_mn (eigenvalue *eval, P_gas p_g, int m, int n);
void eigenvector_mn (eigenvector *evec, P_gas p_g, int m, int n);
void fill_with_vector (double *G, double *V1, double *V2, P_she p_s, int m, int n, double vector[N_]);


#endif //UNTITLED2_HARMONIC_H
