#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "gas_two.h"
#include "residuals.h"
#include "gnuploting.h"

#define IT_max 1

#define SOKOL 1

void print_norms(FILE *f, double *norm, int it_max);

int deg2inx (int x);

int main() {

    P_gas p_d;
    P_she p_s;
    double *G;
    double *V1;
    double *V2;
    int *st;
    int it_sp, it_t, k;
    int it_max = IT_max + 1;
    double *nc_g, *nl2_g, *nw12_g, *nc_v1, *nl2_v1, *nw12_v1, *nc_v2, *nl2_v2, *nw12_v2;
    int deg = deg2inx(it_max+1);
    Norm_Step n_s;

    FILE *f;

    nc_g = (double *) malloc(deg * sizeof(double));
    nl2_g = (double *) malloc(deg * sizeof(double));
    nc_v1 = (double *) malloc(deg * sizeof(double));
    nl2_v1 = (double *) malloc(deg * sizeof(double));
    nc_v2 = (double *) malloc(deg * sizeof(double));
    nl2_v2 = (double *) malloc(deg * sizeof(double));
    nw12_g = (double *) malloc(deg * sizeof(double));
    nw12_v1 = (double *) malloc(deg * sizeof(double));
    nw12_v2 = (double *) malloc(deg * sizeof(double));

    param_dif(&p_d);

    if (SMOOTH_SOLUTION != 1) {
        clean_png("g");
        clean_png("u1");
    }

    k = 0;
    for (it_t = 1; it_t < deg2inx(it_max) + 1; it_t *= 2)
    {
        for (it_sp = 1; it_sp < deg2inx(it_max) + 1; it_sp *= 2)
        {
            param_she_step(&p_s, p_d, it_t, it_sp);

            st = (int*) malloc(p_s.Dim * sizeof(int));
            G = (double*) malloc(p_s.Dim * sizeof(double));
            V1 = (double*) malloc(p_s.Dim * sizeof(double));
            V2 = (double*) malloc(p_s.Dim * sizeof(double));

            Setka(st, &p_s);

            Sxema(G, V1, V2, st, p_s, p_d);

            printf("asdfsdfa \n M_x = %d \n M_y = %d \n N   = %d \n", p_s.M_x, p_s.M_y, p_s.N);

            if (SMOOTH_SOLUTION == 1) {
                residual_Ch(V1, V2, G, p_s, &n_s, u1, u2, g);
                nc_v1[k] = n_s.V1norm;
                nc_v2[k] = n_s.V2norm;
                nc_g[k] = n_s.Gnorm;

                residual_L2h(V1, V2, G, p_s, &n_s, u1, u2, g);
                nl2_v1[k] = n_s.V1norm;
                nl2_v2[k] = n_s.V2norm;
                nl2_g[k] = n_s.Gnorm;

                residual_W12(V1, V2, G, p_s, &n_s, u1, u2, g);
                nw12_v1[k] = n_s.V1norm;
                nw12_v2[k] = n_s.V2norm;
                nw12_g[k] = n_s.Gnorm;
            }
            free(V1);
            free(V2);
            free(G);
//sleep(10);
            ++k;
        }
    }
    if (SMOOTH_SOLUTION == 1) {
        f = fopen("ff.txt", "w");

        if (!f) printf("dsfasafdasdfasdfsafdkajdsfsa \n");

        print_norms(f, nc_g, it_max);
        print_norms(f, nl2_g, it_max);
        print_norms(f, nw12_g, it_max);
        print_norms(f, nc_v1, it_max);
        print_norms(f, nl2_v1, it_max);
        print_norms(f, nw12_v1, it_max);
        print_norms(f, nc_v2, it_max);
        print_norms(f, nl2_v2, it_max);
        print_norms(f, nw12_v2, it_max);

        fclose(f);
    }
    free(nc_g);
    free(nc_v1);
    free(nc_v2);
    free(nl2_g);
    free(nl2_v1);
    free(nl2_v2);
    free(nw12_g);
    free(nw12_v1);
    free(nw12_v2);

    return 0;
}

int deg2inx (int x)
{
    int i;
    int ans = 1;

    if (x < 0) return 0;
    if (x == 1) return 2;
    for (i = 1; i < x+1; ++i)
        ans *= 2;
    return ans;
}

void print_norms(FILE *f, double *norm, int it_max)
{
    int k, it_t, it_sp;

    printf("start print norm \n");
    fprintf(f, "\n\n");
    fprintf(f, "\\begin{tabular}{c r r r r}\n");
    fprintf(f, "\\hline \n");
    fprintf(f, "N \\texttt{\\char`\\\\} M ");
    for (it_sp = 1; it_sp < deg2inx(it_max) + 1; it_sp *= 2)
        fprintf(f, "& %d", zero_spl_x * it_sp);
    fprintf(f, "\\\\ \n\\hline \n");

    k = 0;
    for (it_t = 1; it_t < deg2inx(it_max) + 1; it_t *= 2)
    {
        fprintf(f,"%d ", zero_spl_t * it_t);
        for (it_sp = 1; it_sp < deg2inx(it_max) + 1; it_sp *= 2)
        {
            fprintf (f, "& %e", norm[k]);
            ++k;
        }
        fprintf(f, "\\\\ \n");
    }
    fprintf(f, "\\hline \n\\end{tabular}");
}