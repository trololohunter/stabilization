//
// Created by vover on 4/3/18.
//

#ifndef UNTITLED2_FUNCTIONS_H
#define UNTITLED2_FUNCTIONS_H

typedef double (*func)(double, double, double);

double u1_zero (double t, double x1, double x2);
double u2_zero (double t, double x1, double x2);
double ro_zero (double t, double x1, double x2);
double u1_test_pop (double t, double x1, double x2);
double u2_test_pop (double t, double x1, double x2);
double ro_test_pop (double t, double x1, double x2);
double  g_test_pop (double t, double x1, double x2);
double rhs_f0_test_pop (double t, double x1, double x2);
double rhs_f1_test_pop (double t, double x1, double x2, double mu, double p_ro);
double rhs_f2_test_pop (double t, double x1, double x2, double mu, double p_ro);
double Func_2 (double t, double x, double y, double p_ro, double mu);
double Func_1 (double t, double x, double y, double p_ro, double mu);
double Func_0 (double t, double x, double y);
double g (double t, double x, double y);
double u1 (double t, double x, double y);
double u2 (double t, double x, double y);

#endif //UNTITLED2_FUNCTIONS_H
