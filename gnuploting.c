//
// Created by vover on 6/1/18.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <zconf.h>
#include <stdlib.h>
#include "gnuploting.h"
#include "gas_two.h"

#define LEN 256

#define FNAME_COM "commands.txt"
#define FNAME_G "g.txt"
#define FNAME_U "u.txt"

void print_values_exp (P_she p_s, double *G, FILE *fp);
double get_x_delta (double x, double y, double h);
double get_y_delta (double x, double y, double h);
void print_values (P_she p_s, double *V1, double *V2, FILE *fp);

void run_gnuplot_ (void)
{
    const char *filename_com = FNAME_COM;
    char path_dest [LEN];

    sprintf (path_dest, "\"gnuplot\" %s%c\"", filename_com,'\0');

    int res = system (path_dest);
    (void) res;

    usleep(10);
    return;
}

void fill_files_with_values (P_she p_s, double *V1, double *V2, double *G);
void fill_files_with_values (P_she p_s, double *V1, double *V2, double *G)
{
    FILE *fp_g;
    FILE *fp_u;
    const char *filename_g = FNAME_G;
    const char *filename_u = FNAME_U;

    fp_g = fopen(filename_g, "w");
    fp_u = fopen(filename_u, "w");

    print_values_exp(p_s, G, fp_g);
    print_values(p_s, V1, V2, fp_u);

    fclose(fp_g);
    fclose(fp_u);

    return;
}

void fill_command_file (const char *val_name, P_she p_s, int time_step,
                        const char *val_path)
{
    const char *filename_com = FNAME_COM;
    double t = time_step * p_s.tau;

    FILE *f_com;
    f_com = fopen (filename_com, "w");

    fprintf (f_com, "set terminal png size 1024, 1024\n");
    fprintf (f_com, "set output '%s/%d.png'\n", val_name, time_step);
    fprintf (f_com, "set xlabel \"t = %.4f\" font \"Times-Roman,30\"\n", t);
    fprintf (f_com, "set xrange [0:2*pi]; set yrange [0:2*pi]\n");
    fprintf (f_com, "%s\n", val_path);

    fclose (f_com);
}

void print_paint_pm3d_command (char *path, const char *fname)
{
    const char *palette = "set palette rgbformulae 34,35,36";
    const char *map = "set pm3d map";
    sprintf (path, "%s; %s; splot '%s' with pm3d%c",
             palette, map, fname, '\0');
}

void print_paint_vectors_command (char *path, const char *fname)
{
    const char *palette = "set palette rgbformulae 18,20,22";
    sprintf (path, "plot '%s' using 1:2:3:4 w vec lw 3 filled head%c",
 //   sprintf (path, "%s; plot '%s' using 1:2:($3/(sqrt(($3-$1)**2+($4-$2)**2))):($4/(sqrt(($3-$1)**2+($4-$2)**2))) w vec lw 3 filled head%c",
             fname, '\0');
}

void print_paint_command (char *path, const char *fname)
{
    sprintf (path, "splot '%s'%c", fname, '\0');
}

void run_gnuplot(P_she p_s, double *V1, double *V2, double *G, int time_step)
{
    double t = time_step * p_s.tau;

    fill_files_with_values (p_s, V1, V2, G);
    char path_g [LEN];
    char path_u [LEN];
    print_paint_pm3d_command (path_g, FNAME_G);
    print_paint_vectors_command (path_u, FNAME_U);

    fill_command_file ("g",  p_s, time_step, path_g);
    run_gnuplot_ ();
    fill_command_file ("u1", p_s, time_step, path_u);
    run_gnuplot_ ();

    return;
}

void print_values_exp (P_she p_s, double *G, FILE *fp)
{
    double x_coord = 0;
    double y_coord = 0;
    double hx = p_s.h_x;
    double hy = p_s.h_y;

    for (int i = 0; i < p_s.M_x+1; ++i, x_coord += hx)
    {
        y_coord = 0;
        for (int j = 0; j < p_s.M_y+1; ++j, y_coord += hy)
        {
            fprintf (fp, "%e %e %e\n", x_coord, y_coord, exp (G[j*(p_s.M_x + 1) + i]));
        }
        fprintf (fp, "\n");
    }

}

double get_x_delta (double x, double y, double h)
{
    double d = sqrt (x * x + y * y);
    if (d < 1e-16)
        return 0;
    d /= h;
    return x / d;
}

double get_y_delta (double x, double y, double h)
{
    double d = sqrt (x * x + y * y);
    if (d < 1e-16)
        return 0;
    d /= h;
    return y / d;
}

void print_values (P_she p_s, double *V1, double *V2, FILE *fp)
{
    double x_coord = 0;
    double y_coord = 0;
    double h = p_s.h_x;
    double v1, v2;

    // x y x_delta y_delta
    for (int i = 0; i < p_s.Dim; i++)
    {
        v1 = V1[i];
        v2 = V2[i];
        x_coord = i % (p_s.M_x + 1) * p_s.h_x;
        y_coord = (i / (p_s.M_x + 1)) * p_s.h_y;
        fprintf (fp, "%e %e %e %e\n", x_coord, y_coord,
                 get_x_delta (v1, v2, h), get_y_delta (v1, v2, h));
    }
}

void clean_png (const char *path)
{
    char path_dest [LEN];

    sprintf (path_dest, "rm -rf %s\\/*.png%c\"", path,'\0');

    int res = system (path_dest);
    (void) res;

    usleep(100);
}