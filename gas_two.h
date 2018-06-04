//
// Created by vover on 3/26/18.
//

#ifndef UNTITLED2_GAS_TWO_H
#define UNTITLED2_GAS_TWO_H

#include <math.h>
#include <glob.h>

#define zero_spl_x  20 //zero splitting x
#define zero_spl_y  20
#define zero_spl_t  20

#define SMOOTH_SOLUTION 1
#define SOKOL 1
#define GAMMA 1.4

typedef struct
{

    double Segm_T;
    double Segm_X;
    double Segm_Y;
    double p_ro;
    double mu;
    double omega;

} P_gas;

typedef struct
{

    int M_x;
    int M_y;
    int N;
    int Dim;
    double h_x;
    double h_y;
    double tau;
    double eta;

} P_she;

typedef struct
{

    double thx;
    double thy;
    double thx05;
    double thy05;
    double thx2;
    double thy2;
    double thx4;
    double thy4;
    double thx32;
    double thy32;
    double tau2;
    double tau4;
    double tau6;
    double thxx8;
    double thxx6;
    double thyy8;
    double thyy6;
    double thxy;
    double Max;
    double May;

} T_const;

typedef struct
{

    double MUM;
    double MU8x;
    double MU8y;
    double MU6x;
    double MU6y;
    double MUv1;
    double MUv2;

} MUM_const;

typedef struct
{

    size_t mmg00;
    size_t mmv100;
    size_t mmv200;
    size_t mmgL0;
    size_t mmv1L0;
    size_t mmv2L0;
    size_t mmgR0;
    size_t mmv1R0;
    size_t mmv2R0;
    size_t mmg0L;
    size_t mmv10L;
    size_t mmv20L;
    size_t mmg0R;
    size_t mmv10R;
    size_t mmv20R;

} MM_step;

void param_dif (P_gas *p_d);
void param_she_step(P_she *p_s, P_gas p_d, int it_t, int it_sp);
void param_t_const (T_const *t_c, P_she p_s, P_gas p_d);
void param_MUM_const (MUM_const *MUM_c, P_she p_s, double GG, P_gas p_d);
void param_MM_step (MM_step *MM_s, size_t mm, int n, int m);
void Setka (int *st, P_she *p_s);
void Sxema (double *G, double *V1, double *V2, int *st, P_she p_s, P_gas p_d);
void Sxema_Sokolov (double *G, double *V1, double *V2, int *st, P_she p_s, P_gas p_d);

#endif //UNTITLED2_GAS_TWO_H
