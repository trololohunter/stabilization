//
// Created by vover on 3/22/20.
//

void udxdx (double *u_out, double *u_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y) {
    int i = 0, j = 0;

    for (i = 0; i < N_y; ++i) {
        u_out[i*(N_x+1)+0] = 0;
        u_out[i*(N_x+1)+N_x] = 0;
    }

    for (i = 0; i < N_y; ++i)
        for (j = 1; j < N_x; ++j)
            u_out[i*(N_x+1)+j] = (u_in[i*(N_x+1)+j+1]-2*u_in[i*(N_x+1)+j]+u_in[i*(N_x+1)+j-1])/(h_x*h_x);

    return;
}

void udydy (double *u_out, double *u_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y) {
    int i = 0, j = 0;

    for (j = 0; j < N_x+1; ++j) {
        u_out[j] = (u_in[(N_x+1)+j]-u_in[j])/(h_y*h_y);
        u_out[(N_y-1)*(N_x+1)+j] = (u_in[(N_y-1)*(N_x+1)+j]-u_in[(N_y-2)*(N_x+1)+j])/(h_y*h_y);
    }

    for (i = 1; i < N_y-1; ++i)
        for (j = 0; j < N_x+1; ++j)
            u_out[i*(N_x+1)+j] = (u_in[(i+1)*(N_x+1)+j]-2*u_in[i*(N_x+1)+j]+u_in[(i-1)*(N_x+1)+j])/(h_y*h_y);

    return;
}



void udxdy (double *u_out, double *u_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y) {
    int i = 0, j = 0;

    for (i = 0; i < N_y; ++i) {
        u_out[i*(N_x+1)+0] = 0;
        u_out[i*(N_x+1)+N_x] = 0;
    }

    for (j = 0; j < N_x+1; ++j) {
        u_out[j] = 0;
        u_out[(N_y-1)*(N_x+1)+j] = 0;
    }

    for (i = 1; i < N_y-1; ++i)
        for (j = 1; j < N_x; ++j)
            u_out[i*(N_x+1)+j] = (u_in[(i+1)*(N_x+1)+j+1]-u_in[(i+1)*(N_x+1)+j-1]-u_in[(i-1)*(N_x+1)+j+1]+u_in[(i-1)*(N_x+1)+j-1])/(h_x*h_y);

    return;
}

void vdxdx (double *v_out, double *v_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y) {
    int i = 0, j = 0;

    for (i = 0; i < N_y+1; ++i) {
        v_out[i*N_x + 0] = (v_in[i*N_x + 0]-v_in[i*N_x + 1])/(h_x*h_x);
        v_out[i*N_x+N_x-1] = (v_in[i*N_x+N_x-2]-v_in[i*N_x+N_x-1])/(h_x*h_x);
    }

    for (i = 0; i < N_y+1; ++i)
        for (j = 1; j < N_x-1; ++j)
            v_out[i*(N_x)+j] = (v_in[i*(N_x)+j+1]-2*v_in[i*(N_x)+j]+v_in[i*(N_x)+j-1])/(h_x*h_x);

    return;
}

void vdydy (double *v_out, double *v_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y) {
    int i = 0, j = 0;

    for (j = 0; j < N_x; ++j) {
        v_out[j] = 0;
        v_out[N_y*N_x+j] = 0;
    }

    for (i = 1; i < N_y; ++i)
        for (j = 0; j < N_x; ++j)
            v_out[i*N_x+j] = (v_in[(i+1)*N_x+j]-2*v_in[i*N_x+j]+v_in[(i-1)*N_x+j])/(h_y*h_y);

    return;
}



void vdxdy (double *v_out, double *v_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y) {
    int i = 0, j = 0;

    for (i = 0; i < N_y+1; ++i) {
        v_out[i*N_x+0] = 0;
        v_out[i*N_x+N_x-1] = 0;
    }

    for (j = 0; j < N_x; ++j) {
        v_out[j] = 0;
        v_out[N_y*N_x+j] = 0;
    }

    for (i = 1; i < N_y; ++i)
        for (j = 1; j < N_x-1; ++j)
            v_out[i*N_x+j] = (v_in[(i+1)*N_x+j+1]-v_in[(i+1)*N_x+j-1]-v_in[(i-1)*N_x+j+1]+v_in[(i-1)*N_x+j-1])/(h_x*h_y);

    return;
}

void udx (double *p_out, double *u_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y) {
    int i = 0, j = 0;

    for (i = 0; i < N_y; ++i)
        for (j = 0; j < N_x; ++j)
            p_out[i*N_x+j] = (u_in[i*(N_x+1)+j+1]-u_in[i*(N_x+1)+j])/h_x;

    return;
}

void vdy (double *p_out, double *v_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y){
    int i = 0, j = 0;

    for (i = 0; i < N_y; ++i)
        for (j = 0; j < N_x; ++j)
            p_out[i*N_x+j] = (v_in[(i+1)*N_x+j]-v_in[i*N_x+j])/h_y;

    return;
}

void pdx (double *u_out, double *p_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y){
    int i = 0, j = 0;

    for (i = 0; i < N_y; ++i) {
        u_out[i*(N_x+1)+0] = 0;
        u_out[i*(N_x+1)+N_x] = 0;
    }
    for (i = 0; i < N_y; ++i)
        for (j = 1; j < N_x; ++j)
            u_out[i*(N_x+1)+j] = (p_in[i*N_x+j]-p_in[i*N_x+j-1])/h_x;

    return;
}

void pdy (double *v_out, double *p_in, int N_x, int N_y, double h_x, double h_y, double L_x, double L_y){
    int i = 0, j = 0;

    for (j = 0; j < N_x; ++j) {
        v_out[j] = 0;
        v_out[N_y*N_x+j] = 0;
    }
    for (i = 1; i < N_y; ++i)
        for (j = 0; j < N_x; ++j)
            v_out[i*N_x+j] = (p_in[i*N_x+j]-p_in[(i-1)*N_x+j])/h_y;

    return;
}