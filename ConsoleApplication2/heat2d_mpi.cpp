#define MSMPI_NO_SAL
#define _Out_writes_(x)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


static double bc_left = 10.0;   /* u(t, 0,   y) */
static double bc_right = 40.0;   /* u(t, 1,   y) */
static double bc_front = 30.0;   /* u(t, x,   0) */
static double bc_back = 50.0;   /* u(t, x,   1) */

static inline double f_init(int /*i*/, int /*j*/) { return 0.0; }

static inline double* alloc2d(int rows, int cols)
{
    double* p = (double*)calloc((size_t)rows * cols, sizeof(double));
    if (!p) { perror("calloc"); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }
    return p;
}

#define IDX(i, j, cols) ((size_t)(i) * (cols) + (j))

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int myid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int n = (argc > 1) ? atoi(argv[1]) : 62;
    int iterations = (argc > 2) ? atoi(argv[2]) : 10000;
    if (n < 1) n = 62;
    if (iterations < 1) iterations = 10000;

    const int N = n + 2; 
    const double c = 0.1;           
    const double ds = 1.0 / (n + 1);
    const double dt = (ds * ds) / (4.0 * c);
    const double r = c * dt / (ds * ds);

    int base_rows = N / nprocs;
    int extra = N % nprocs;

    int local_rows = base_rows + (myid < extra ? 1 : 0);
    int global_row0 = myid * base_rows + (myid < extra ? myid : extra);

    int lr = local_rows;
    int cols = N;
    int local_total_rows = lr + 2;

    double* u_old = alloc2d(local_total_rows, cols);
    double* u_new = alloc2d(local_total_rows, cols);


    for (int li = 0; li <= lr + 1; li++) {
        int gi = global_row0 - 1 + li;
        for (int j = 0; j < cols; j++) {
            double val;
            if (gi == 0)     val = bc_back;
            else if (gi == N - 1) val = bc_front;
            else if (j == 0)      val = bc_left;
            else if (j == N - 1)  val = bc_right;
            else                  val = f_init(gi, j);
            u_old[IDX(li, j, cols)] = val;
            u_new[IDX(li, j, cols)] = val;
        }
    }

    int ifirst = 1;
    int ilast = lr;
    if (myid == 0)          ifirst++;
    if (myid == nprocs - 1) ilast--;

    MPI_Datatype row_type;
    MPI_Type_contiguous(cols, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();


    for (int iter = 0; iter < iterations; iter++) {

        for (int li = ifirst; li <= ilast; li++) {
            for (int j = 1; j < cols - 1; j++) {
                u_new[IDX(li, j, cols)] =
                    u_old[IDX(li, j, cols)]
                    + r * (u_old[IDX(li + 1, j, cols)]
                        + u_old[IDX(li - 1, j, cols)]
                        - 4.0 * u_old[IDX(li, j, cols)]
                        + u_old[IDX(li, j + 1, cols)]
                        + u_old[IDX(li, j - 1, cols)]);
            }
        }

        MPI_Request reqs[4];
        int nreqs = 0;

        if (myid < nprocs - 1) {
            MPI_Isend(&u_new[IDX(lr, 0, cols)], 1, row_type,
                myid + 1, 0, MPI_COMM_WORLD, &reqs[nreqs++]);
            MPI_Irecv(&u_new[IDX(lr + 1, 0, cols)], 1, row_type,
                myid + 1, 1, MPI_COMM_WORLD, &reqs[nreqs++]);
        }
        if (myid > 0) {
            MPI_Isend(&u_new[IDX(1, 0, cols)], 1, row_type,
                myid - 1, 1, MPI_COMM_WORLD, &reqs[nreqs++]);
            MPI_Irecv(&u_new[IDX(0, 0, cols)], 1, row_type,
                myid - 1, 0, MPI_COMM_WORLD, &reqs[nreqs++]);
        }
        if (nreqs > 0)
            MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);

        /* ── 3. swap buffers ── */
        double* tmp = u_old;
        u_old = u_new;
        u_new = tmp;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed = MPI_Wtime() - t_start;

    int* rcounts = NULL, * displs = NULL;
    double* full_grid = NULL;

    if (myid == 0) {
        rcounts = (int*)malloc(nprocs * sizeof(int));
        displs = (int*)malloc(nprocs * sizeof(int));
        full_grid = alloc2d(N, cols)

        int offset = 0;
        for (int p = 0; p < nprocs; p++) {
            int pr = base_rows + (p < extra ? 1 : 0);

            int owned;
            if (p == 0 && nprocs == 1) owned = pr - 2;
            else if (p == 0 || p == nprocs - 1) owned = pr - 1;
            else                                owned = pr;
            rcounts[p] = owned * cols;
            displs[p] = offset * cols;
            offset += owned;
        }
    }

    int owned_rows;
    if (nprocs == 1)          owned_rows = lr - 2;
    else if (myid == 0 || myid == nprocs - 1) owned_rows = lr - 1;
    else                           owned_rows = lr;

    int send_start = ifirst;
    double* send_buf = &u_old[IDX(send_start, 0, cols)];

    MPI_Gatherv(send_buf, owned_rows * cols, MPI_DOUBLE,
        full_grid, rcounts, displs, MPI_DOUBLE,
        0, MPI_COMM_WORLD);

    if (myid == 0) {
        printf("========================================\n");
        printf("  2D Heat Equation — Parallel MPI Solver\n");
        printf("  (Horak & Gruber, Parallel Numerics '05)\n");
        printf("========================================\n");
        printf("  Grid:          %d × %d  (n = %d inner pts)\n", N, N, n);
        printf("  Processors:    %d\n", nprocs);
        printf("  Iterations:    %d\n", iterations);
        printf("  c (diffusivity): %.4f\n", c);
        printf("  ds:            %.6f\n", ds);
        printf("  dt:            %.6e\n", dt);
        printf("  Fourier no. r: %.6f  (stable if <= 0.25)\n", r);
        printf("  Wall time:     %.6f s\n", elapsed);
        printf("----------------------------------------\n");

        printf("  Sample temperatures (inner grid):\n");
        int pts[][2] = { {N / 4, N / 4}, {N / 2, N / 2}, {3 * N / 4, 3 * N / 4},
                        {N / 4, 3 * N / 4}, {3 * N / 4, N / 4} };
        int npts = (int)(sizeof(pts) / sizeof(pts[0]));
        for (int k = 0; k < npts; k++) {
            int gi = pts[k][0], gj = pts[k][1];
            int fi = gi - 1;
            if (fi >= 0 && fi < N - 2)
                printf("    u(x=%d/%d, y=%d/%d) = %.4f\n",
                    gi, N - 1, gj, N - 1,
                    full_grid[IDX(fi, gj, cols)]);
        }

        double umin = 1e30, umax = -1e30;
        for (int i = 0; i < N - 2; i++) {
            for (int j = 1; j < cols - 1; j++) {
                double v = full_grid[IDX(i, j, cols)];
                if (v < umin) umin = v;
                if (v > umax) umax = v;
            }
        }
        printf("  Interior  min = %.4f,  max = %.4f\n", umin, umax);
        printf("========================================\n");

        free(rcounts); free(displs); free(full_grid);
    }

    double max_time;
    MPI_Reduce(&elapsed, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myid == 0)
        printf("  Max wall time across all ranks: %.6f s\n\n", max_time);

    MPI_Type_free(&row_type);
    free(u_old); free(u_new);

    MPI_Finalize();
    return 0;
}