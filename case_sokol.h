//
// Created by vover on 6/4/18.
//

#ifndef UNTITLED2_CASE_SOKOL_H
#define UNTITLED2_CASE_SOKOL_H

#include "functions.h"


void first_fill_sokol (double *V1, double *V2, double *G, P_she p_s,
                 double w, func u1, func u2, func ro);
void first_fill_sokol___ (double *V1, double *V2, double *G, P_she p_s,
                    double w);
void case_sokol_H (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                   double *V1, double *V2, double *G,
                   int m, P_she p_s, double mu, size_t mm);
void case0_sokol_V (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm, double p_ro);
void case1_sokol_V (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm);
void case2_sokol_V (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm);
void case3_sokol_V (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm, double w);
void case4_sokol_V (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm);
void case9_sokol_V (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                    double *V1, double *V2, double *G, double *Hnew,
                    int m, P_she p_s, double mu, size_t mm);
void case10_sokol_V (QMatrix_L *H, Vector *b_H, T_const t_c, int k,
                     double *V1, double *V2, double *G, double *Hnew,
                     int m, P_she p_s, double mu, size_t mm);

#endif //UNTITLED2_CASE_SOKOL_H
