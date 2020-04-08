//
// Created by vover on 4/3/18.
//

#include <math.h>
#include "functions.h"
#include "gas_two.h"

double u1_zero (double t, double x1, double x2)
{
    (void) t;
    (void) x1;
    (void) x2;
    return 0;
}

double u2_zero (double t, double x1, double x2)
{
    (void) t;
    (void) x1;
    (void) x2;
    return 0;
}

double ro_zero (double t, double x1, double x2)
{
    (void) t;
    (void) x1;
    (void) x2;
    return 1;
}

double u1_test_pop (double t, double x1, double x2)
{
    //return 0;
    return sin(2 * M_PI * x1) * sin(2 * M_PI * x2) * exp (t);
}

double u2_test_pop (double t, double x1, double x2)
{
    //return 0;
    return sin(2 * M_PI * x1) * sin(2 * M_PI * x2) * exp (-t);
}

double ro_test_pop (double t, double x1, double x2)
{
    //return 1;
    return (cos(2 * M_PI * x1) + 1.5) * (sin(2 * M_PI * x2) + 1.5) * exp(t);
}

double  g_test_pop (double t, double x1, double x2)
{
    return log(ro_test_pop(t, x1, x2));
}
/*
double du1dt (double t, double x1, double x2)
{
    return u1_test_pop(t, x1, x2);
}

double du2dt (double t, double x1, double x2)
{
    return -u2_test_pop(t, x1, x2);
}

double dgdt (double t, double x1, double x2)
{
    return 1;
}

double du1dx1 (double t, double x1, double x2)
{
    return 2 * M_PI * exp(t) * cos (2 * M_PI * x1) * sin (2 * M_PI * x2);
}

double du2dx1 (double t, double x1, double x2)
{
    return 2 * M_PI * exp(-t) * cos (2 * M_PI * x1) * sin (2 * M_PI * x2);
}

double dgdx1 (double t, double x1, double x2)
{
    return -2 * M_PI * sin (2 * M_PI * x1) / (cos(2 * M_PI * x1) + 1.5);
}

double du1dx2 (double t, double x1, double x2)
{
    return 2 * M_PI * exp(t) * sin (2 * M_PI * x1) * cos (2 * M_PI * x2);
}

double du2dx2 (double t, double x1, double x2)
{
    return 2 * M_PI * exp(-t) * sin (2 * M_PI * x1) * cos (2 * M_PI * x2);
}

double dgdx2 (double t, double x1, double x2)
{
    return 2 * M_PI * cos(2 * M_PI * x2)/(sin(2 * M_PI * x2) + 1.5);
}

double du1u2dx1 (double t, double x1, double x2)
{
    return 2 * M_PI * sin (4 * M_PI * x1) * sin (2 * M_PI * x2) * sin (2 * M_PI * x2);
}

double du1u2dx2 (double t, double x1, double x2)
{
    return 2 * M_PI * sin (4 * M_PI * x2) * sin (2 * M_PI * x1) * sin (2 * M_PI * x1);
}

double du1u1dx1 (double t, double x1, double x2)
{
    return 2 * M_PI * sin (4 * M_PI * x1) * sin (2 * M_PI * x2) * sin (2 * M_PI * x2) * exp(2 * t);
}

double du2u2dx2 (double t, double x1, double x2)
{
    return 2 * M_PI * sin (4 * M_PI * x2) * sin (2 * M_PI * x1) * sin (2 * M_PI * x1) * exp(-2 * t);
}

double du1gdx1 (double t, double x1, double x2)
{
    return 2 * M_PI * exp(t) * cos(2 * M_PI * x1) * sin(2 * M_PI * x2) * g_test_pop(t, x1, x2)
           - 2 * M_PI * exp(t) * sin (2 * M_PI * x1) * sin (2 * M_PI * x1) * sin (2 * M_PI * x2)/(cos(2 * M_PI * x1) + 1.5);
}

double du2gdx2 (double t, double x1, double x2)
{
    return 2 * M_PI * exp(-t) * sin(2 * M_PI * x1) * cos(2 * M_PI * x2) * g_test_pop(t, x1, x2)
           + 2 * M_PI * exp(-t) * sin(2 * M_PI * x1) * sin (2 * M_PI * x2) * cos (2 * M_PI * x2) / (sin(2 * M_PI * x2) + 1.5);
}

double d2u1dx1dx1 (double t, double x1, double x2)
{
    return -4 * M_PI * M_PI * exp(t) * sin(2 * M_PI * x1) * sin(2 * M_PI * x2);
}

double d2u1dx1dx2 (double t, double x1, double x2)
{
    return 4 * M_PI * M_PI * exp(t) * cos(2 * M_PI * x1) * cos(2 * M_PI * x2);
}

double d2u1dx2dx2 (double t, double x1, double x2)
{
    return -4 * M_PI * M_PI * exp(t) * sin(2 * M_PI * x1) * sin(2 * M_PI * x2);
}

double d2u2dx2dx2 (double t, double x1, double x2)
{
    return -4 * M_PI * M_PI * exp(-t) * sin(2 * M_PI * x1) * sin(2 * M_PI * x2);
}

double d2u2dx1dx2 (double t, double x1, double x2)
{
    return 4 * M_PI * M_PI * exp(-t) * cos(2 * M_PI * x1) * cos(2 * M_PI * x2);
}

double d2u2dx1dx1 (double t, double x1, double x2)
{
    return -4 * M_PI * M_PI * exp(-t) * sin(2 * M_PI * x1) * sin(2 * M_PI * x2);
}

double rhs_f0_test_pop (double t, double x1, double x2)
{
//    return 0;
    return dgdt(t, x1, x2)
           + 1./2. * (u1_test_pop(t,x1,x2) * dgdx1(t,x1,x2)
                      + du1gdx1(t,x1,x2)
                      + (2 - g_test_pop(t,x1,x2)) * du1dx1(t,x1,x2))
           + 1./2. * (u2_test_pop(t,x1,x2) * dgdx2(t,x1,x2)
                      + du2gdx2(t,x1,x2)
                      + (2 - g_test_pop(t,x1,x2)) * du2dx2(t,x1,x2));
}

double rhs_f1_test_pop (double t, double x1, double x2, double mu, double p_ro)
{//return 0;
    return du1dt(t,x1,x2)
            + 1./3. * (u1_test_pop(t,x1,x2) * du1dx1(t,x1,x2) + du1u1dx1(t,x1,x2))
            + 1./2. * (u2_test_pop(t,x1,x2) * du1dx2(t,x1,x2) + du1u2dx2(t,x1,x2)
                       - u1_test_pop(t,x1,x2) * du2dx2(t,x1,x2))
            + p_ro * dgdx1(t,x1,x2)
            - mu/ro_test_pop(t,x1,x2)
            * (4./3. * d2u1dx1dx1(t,x1,x2) + d2u1dx2dx2(t,x1,x2)
               + 1./3. * d2u2dx1dx2(t,x1,x2));

}

double rhs_f2_test_pop (double t, double x1, double x2, double mu, double p_ro)
{//return 0;

    return du2dt(t,x1,x2)
            + 1./3. * (u2_test_pop(t,x1,x2) * du2dx2(t,x1,x2) + du2u2dx2(t,x1,x2))
            + 1./2. * (u1_test_pop(t,x1,x2) * du2dx1(t,x1,x2) + du1u2dx1(t,x1,x2)
                       - u2_test_pop(t,x1,x2) * du1dx1(t,x1,x2))
            + p_ro * dgdx2 (t,x1,x2)
            - mu/ro_test_pop(t,x1,x2)
            * (4./3. * d2u2dx2dx2(t,x1,x2) + d2u2dx1dx1(t,x1,x2)
               + 1./3. * d2u1dx1dx2(t,x1,x2));

}
*/

double u (double coeff, double x, double y, int m, int n)
{
    return sin(m*x/2) * cos(n*y/2) * coeff;
}

double v (double coeff, double x, double y, int m, int n)
{
    return cos(m*x/2) * sin(n*y/2) * coeff;
}

double p (double coeff, double x, double y, int m, int n)
{
    return cos(m*x/2) * cos(n*y/2) * coeff;
}

double u1 (double t, double x, double y)
{
    //double res = sin (x) * sin (y) * exp (t);
    //return res;
    return 0;
}

double u2 (double t, double x, double y)
{
    //double res = sin (x) * sin (y) * exp (-t);
    //return res;
    return 0;
}

double rho (double t, double x, double y)
{
    //printf ("Do not use rho (t, x, y)!\n");
    double res = (cos (x) + 3./2.)
                 * (sin (y) + 3./2.) * exp (t);
    return res;
}

double g (double t, double x, double y)
{
    return 1;
    double res = (cos (x) + 3./2.) * (sin (y) + 3./2.);
    res = log (res);
    res += t;
    return res;
}

///////////////////////////////////////////////////////////////////////////////

double dg_dt (double t, double x, double y)
{
    (void) t;
    (void) x;
    (void) y;

    return 1;
}

double dg_dx (double t, double x, double y)
{
    (void) t;
    (void) y;

    double res  = - sin (x);
    res /= (cos ( x) + 3./2.);
    return res;
}

double dg_dy (double t, double x, double y)
{
    (void) t;
    (void) x;

    double res = cos (y);
    res /=  (sin (y) + 3./2.);
    return res;
}

double du1_dt (double t, double x, double y)
{
    return u1 (t, x, y);
}

double du1_dx (double t, double x, double y)
{

    double res = cos (x) * sin (y) * exp (t);
    return res;
}

double du1_dy (double t, double x, double y)
{
    double res = sin (x) * cos (y) * exp (t);
    return res;
}

double ddu1_dxdx (double t, double x, double y)
{
    return - u1 (t, x, y);
}

double ddu1_dydy (double t, double x, double y)
{
    return - u1 (t, x, y);
}

double ddu1_dxdy (double t, double x, double y)
{
    return cos (x) * cos (y) * exp (t);
}

double du1u1_dx (double t, double x, double y)
{
    return 2 * u1 (t, x, y) * du1_dx (t, x, y);
}

double du2_dt (double t, double x, double y)
{
    return - u2 (t, x, y);
}

double du2_dx (double t, double x, double y)
{
    double res = cos (x) * sin (y) * exp (-t);
    return res;
}

double du2_dy (double t, double x, double y)
{
    double res = sin (x) * cos (y) * exp (-t);
    return res;
}

double ddu2_dxdx (double t, double x, double y)
{
    return - u2 (t, x, y);
}

double ddu2_dydy (double t, double x, double y)
{
    return - u2 (t, x, y);
}

double ddu2_dxdy (double t, double x, double y)
{
    return cos (x) * cos (y) * exp (-t);
}

double du2u2_dy (double t, double x, double y)
{
    return 2 * u2 (t, x, y) * du2_dy (t, x, y);
}

double du1g_dx (double t, double x, double y)
{
    return u1 (t, x, y) * dg_dx (t, x, y) + du1_dx (t, x, y) * g (t, x, y);
}

double du2g_dy (double t, double x, double y)
{
    return u2 (t, x, y) * dg_dy (t, x, y) + du2_dy (t, x, y) * g (t, x, y);
}

double du1u2_dx (double t, double x, double y)
{
    return u1 (t, x, y) * du2_dx (t, x, y) + du1_dx (t, x, y) * u2 (t, x, y);
}

double du1u2_dy (double t, double x, double y)
{
    return u1 (t, x, y) * du2_dy (t, x, y) + du1_dy (t, x, y) * u2 (t, x, y);
}

///////////////////////////////////////////////////////////////////////////////

double Func_0 (double t, double x, double y)
{
    return 0;
    if (SMOOTH_SOLUTION != 1) return 0;
    double tmp = 2 - g (t, x, y);
    double res =
            + dg_dt (t, x, y)
            + 0.5 * (
                    + u1 (t, x, y) * dg_dx (t, x, y)
                    + du1g_dx (t, x, y)
                    + tmp * du1_dx (t, x, y))
            + 0.5 * (
                    + u2 (t, x, y) * dg_dy (t, x, y)
                    + du2g_dy (t, x, y)
                    + tmp * du2_dy (t, x, y));

    return res;
}

double Func_1 (double t, double x, double y, double p_rho, double mu)
{
    return 0;
    if (SMOOTH_SOLUTION != 1) return 0;
    double res =
            + du1_dt (t, x, y)
            + (1. / 3.) * (u1 (t, x, y) * du1_dx (t, x, y)
                           + du1u1_dx (t, x, y))
            + (1. / 2.) * (u2 (t, x, y) * du1_dy (t, x, y)
                           + du1u2_dy (t, x, y)
                           - u1 (t, x, y) * du2_dy (t, x, y))
            + p_rho * dg_dx (t, x, y)
            - (mu / exp (g (t, x, y))) * ((4. / 3.) * ddu1_dxdx (t, x, y)
                                          + ddu1_dydy (t, x, y)
                                          + (1. / 3.) * (ddu2_dxdy (t, x, y)));

    return res;
}

double Func_2 (double t, double x, double y, double p_rho, double mu)
{
    return 0;
    if (SMOOTH_SOLUTION != 1) return 0;
    double res =
            + du2_dt (t, x, y)
            + (1. / 3.) * (u2 (t, x, y) * du2_dy (t, x, y)
                           + du2u2_dy (t, x ,y))
            + (1. / 2.) * (u1 (t, x, y) * du2_dx (t, x, y)
                           + du1u2_dx (t, x, y)
                           - u2 (t, x, y) * du1_dx (t, x, y))
            + p_rho * dg_dy (t, x, y)
            - (mu / exp (g (t, x, y))) * ((4. / 3.) * ddu2_dydy (t, x, y)
                                          + ddu2_dxdx (t, x, y)
                                          + (1. / 3.) * (ddu1_dxdy (t, x, y)));

    return res;
}

