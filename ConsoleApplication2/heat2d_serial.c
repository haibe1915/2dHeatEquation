#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define BC_LEFT   10.0
#define BC_RIGHT  40.0
#define BC_FRONT  30.0
#define BC_BACK   50.0

#define IC_VALUE  0.0

static double *alloc2d(int rows, int cols)
{
    double *p = (double *)calloc((size_t)rows * cols, sizeof(double));
    if (!p) { perror("calloc"); exit(EXIT_FAILURE); }
    return p;
}

#define IDX(i, j, cols) ((size_t)(i) * (cols) + (j))

static double get_time(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char *argv[])
{
    int n          = (argc > 1) ? atoi(argv[1]) : 62;
    int iterations = (argc > 2) ? atoi(argv[2]) : 10000;
    if (n < 1) n = 62;
    if (iterations < 1) iterations = 10000;

    const int    N  = n + 2;
    const double c  = 0.1;
    const double ds = 1.0 / (n + 1);
    const double dt = (ds * ds) / (4.0 * c);
    const double r  = c * dt / (ds * ds);

    double *u_old = alloc2d(N, N);
    double *u_new = alloc2d(N, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double val;

            if      (i == 0)     val = BC_BACK;
            else if (i == N - 1) val = BC_FRONT;
            else if (j == 0)     val = BC_LEFT;
            else if (j == N - 1) val = BC_RIGHT; 
            else                 val = IC_VALUE;  

            u_old[IDX(i, j, N)] = val;
            u_new[IDX(i, j, N)] = val;
        }
    }

    double t_start = get_time();

    for (int iter = 0; iter < iterations; iter++) {

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                u_new[IDX(i, j, N)] =
                    u_old[IDX(i,     j,     N)]
                  + r * ( u_old[IDX(i + 1, j,     N)]
                         + u_old[IDX(i - 1, j,     N)]
                         - 4.0 * u_old[IDX(i, j,   N)]
                         + u_old[IDX(i,     j + 1, N)]
                         + u_old[IDX(i,     j - 1, N)] );
            }
        }

        double *tmp = u_old;
        u_old = u_new;
        u_new = tmp;
    }

    double elapsed = get_time() - t_start;


    printf("========================================\n");
    printf("  2D Heat Equation — Serial C Solver\n");
    printf("  (Horak & Gruber, Parallel Numerics '05)\n");
    printf("========================================\n");
    printf("  Grid:            %d x %d  (n = %d inner pts)\n", N, N, n);
    printf("  Processors:      1  (serial)\n");
    printf("  Iterations:      %d\n", iterations);
    printf("  c (diffusivity): %.4f\n", c);
    printf("  ds:              %.6f\n", ds);
    printf("  dt:              %.6e\n", dt);
    printf("  Fourier no. r:   %.6f  (stable if <= 0.25)\n", r);
    printf("  Wall time:       %.6f s\n", elapsed);
    printf("----------------------------------------\n");

    printf("  Sample temperatures (inner grid):\n");
    int pts[][2] = {
        {N/4,     N/4    },
        {N/2,     N/2    },
        {3*N/4,   3*N/4  },
        {N/4,     3*N/4  },
        {3*N/4,   N/4    }
    };
    int npts = (int)(sizeof(pts) / sizeof(pts[0]));
    for (int k = 0; k < npts; k++) {
        int gi = pts[k][0], gj = pts[k][1];
        if (gi > 0 && gi < N-1 && gj > 0 && gj < N-1)
            printf("    u(x=%d/%d, y=%d/%d) = %.4f\n",
                   gi, N-1, gj, N-1,
                   u_old[IDX(gi, gj, N)]);
    }

    double umin =  1e30;
    double umax = -1e30;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            double v = u_old[IDX(i, j, N)];
            if (v < umin) umin = v;
            if (v > umax) umax = v;
        }
    }
    printf("  Interior  min = %.4f,  max = %.4f\n", umin, umax);
    printf("========================================\n\n");

    free(u_old);
    free(u_new);
    return 0;
}
