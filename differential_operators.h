//
// Created by vover on 3/22/20.
//

#ifndef UNTITLED2_DIFFERENTIAL_OPERATORS_H
#define UNTITLED2_DIFFERENTIAL_OPERATORS_H

void udxdx (double *u_out, double *u_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void udydy (double *u_out, double *u_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void udxdy (double *u_out, double *u_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void vdxdx (double *v_out, double *v_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void vdydy (double *v_out, double *v_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void vdxdy (double *v_out, double *v_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void udx (double *p_out, double *u_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void vdy (double *p_out, double *v_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void pdx (double *u_out, double *p_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);
void pdy (double *v_out, double *p_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y);

#endif //UNTITLED2_DIFFERENTIAL_OPERATORS_H
