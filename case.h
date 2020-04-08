//
// Created by vover on 4/4/18.
//

#ifndef UNTITLED2_CASE_H
#define UNTITLED2_CASE_H

#include "functions.h"

void first_fill (double *V1, double *V2, double *G, P_she p_s,
                 double w, func u1, func u2, func ro);
void first_fill___ (double *V1, double *V2, double *G, P_she p_s,
                 double w);
size_t case0 (QMatrix_L *A, Vector *B, T_const t_c, MUM_const m_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, double p_ro, size_t mm);
size_t case1 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm);
size_t case2 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm);
size_t case3 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm, double w);
size_t case4 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm);
size_t case9 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm);
size_t case10 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *G, int m, P_she p_s, double mu, size_t mm);

#endif //UNTITLED2_CASE_H
