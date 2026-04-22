/*
 * Serial Numerical Solution of 2D Heat Equation
 * Baseline (single-processor) version for comparison with MPI parallel code.
 *
 * Based on:
 *   Horak, V. & Gruber, P. — "Parallel Numerical Solution of 2-D Heat Equation"
 *   Parallel Numerics '05, pp. 47-56
 *
 * Compile:  gcc -O2 -o heat2d_serial heat2d_serial.c -lm
 * Run:      ./heat2d_serial [n] [iterations]
 *
 *   n          = number of inner grid points per side (default 62 → 64x64 grid)
 *   iterations = number of time steps              (default 10000)
 *
 * PDE:      u_t = c * (u_xx + u_yy)
 * Domain:   unit square [0,1]x[0,1]
 * BC:       left=10, right=40, front=30, back=50
 * IC:       all inner points = 0
 * Stability: dt = ds^2 / (4*c)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ── boundary conditions (same as MPI version / paper §3) ── */
#define BC_LEFT   10.0
#define BC_RIGHT  40.0
#define BC_FRONT  30.0
#define BC_BACK   50.0

/* ── initial condition for inner points ── */
#define IC_VALUE  0.0

/* ── helper: allocate flat 2D array (rows x cols) ── */
static double *alloc2d(int rows, int cols)
{
    double *p = (double *)calloc((size_t)rows * cols, sizeof(double));
    if (!p) { perror("calloc"); exit(EXIT_FAILURE); }
    return p;
}

/* ── index macro ── */
#define IDX(i, j, cols) ((size_t)(i) * (cols) + (j))

/* ── wall-clock timer (seconds) ── */
static double get_time(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char *argv[])
{
    /* ── parse arguments ── */
    int n          = (argc > 1) ? atoi(argv[1]) : 62;
    int iterations = (argc > 2) ? atoi(argv[2]) : 10000;
    if (n < 1) n = 62;
    if (iterations < 1) iterations = 10000;

    /*
     * Grid layout (paper §2):
     *   Total points per side N = n + 2  (indices 0 … n+1)
     *   Spatial step  ds = 1 / (n+1)
     *   Time step     dt = ds^2 / (4*c)   ← stability choice from paper §3
     *   Fourier no.   r  = c*dt/ds^2 = 0.25  (stability limit)
     */
    const int    N  = n + 2;
    const double c  = 0.1;
    const double ds = 1.0 / (n + 1);
    const double dt = (ds * ds) / (4.0 * c);
    const double r  = c * dt / (ds * ds);

    /* ── allocate two N×N grids ── */
    double *u_old = alloc2d(N, N);
    double *u_new = alloc2d(N, N);

    /* ──────────────────────────────────────────────
     * Initialise: boundary conditions + initial values
     * (paper §2, single-processor pseudocode)
     * ────────────────────────────────────────────── */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double val;

            if      (i == 0)     val = BC_BACK;    /* top row    */
            else if (i == N - 1) val = BC_FRONT;   /* bottom row */
            else if (j == 0)     val = BC_LEFT;    /* left col   */
            else if (j == N - 1) val = BC_RIGHT;   /* right col  */
            else                 val = IC_VALUE;   /* inner pts  */

            u_old[IDX(i, j, N)] = val;
            u_new[IDX(i, j, N)] = val;
        }
    }

    /* ──────────────────────────────────────────────
     * Main time-stepping loop  (paper §2, serial version)
     *
     *   u_new[i][j] = u_old[i][j]
     *               + r * ( u_old[i+1][j] + u_old[i-1][j]
     *                      - 4*u_old[i][j]
     *                      + u_old[i][j+1] + u_old[i][j-1] )
     *
     * Only inner points (1 ≤ i,j ≤ n) are updated;
     * boundary rows/columns stay fixed.
     * ────────────────────────────────────────────── */
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

        /* swap pointers — no copy needed */
        double *tmp = u_old;
        u_old = u_new;
        u_new = tmp;
    }

    double elapsed = get_time() - t_start;

    /* ──────────────────────────────────────────────
     * Print results
     * ────────────────────────────────────────────── */
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

    /* sample a few interior temperatures */
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

    /* min / max over interior */
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

    /* ──────────────────────────────────────────────
     * Optional: write full grid to CSV for plotting
     * Uncomment if you want to visualise in Excel/Python/MATLAB
     * ────────────────────────────────────────────── */
    /*
    FILE *fp = fopen("heat2d_result.csv", "w");
    if (fp) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(fp, "%.4f", u_old[IDX(i, j, N)]);
                if (j < N-1) fprintf(fp, ",");
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
        printf("  Result saved to heat2d_result.csv\n\n");
    }
    */

    free(u_old);
    free(u_new);
    return 0;
}
