//
// Created by vover on 4/3/18.
//

#include "functions.h"
#include "gas_two.h"

#ifndef UNTITLED2_RESIDUALS_H
#define UNTITLED2_RESIDUALS_H

typedef struct
{
    double V1norm;
    double V2norm;
    double Gnorm;

} Norm_Step;

double residual_Ch_step (double *V1, double *V2, double *G, P_she p_s, int k, func u1, func u2, func ro);
double residual_Ch (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro);
double residual_L2h_step (double *V1, double *V2, double *G, P_she p_s, int k, func u1, func u2, func ro);
double residual_L2h (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro);
double residual_W12 (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro);

double residual_Ch_Sokol (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro);
double residual_L2h_Sokol (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro);
double residual_W12_Sokol (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro);
#endif //UNTITLED2_RESIDUALS_H
