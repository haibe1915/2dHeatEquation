#define MSMPI_NO_SAL
#define _Out_writes_(x)
/*
 * Parallel Numerical Solution of 2D Heat Equation
 * Implementation based on:
 *   Horak, V. & Gruber, P. — "Parallel Numerical Solution of 2-D Heat Equation"
 *   Parallel Numerics '05, pp. 47–56
 *
 * Compile:  mpicc -O2 -o heat2d_mpi heat2d_mpi.c -lm
 * Run:      mpirun -np <num_procs> ./heat2d_mpi [n] [iterations]
 *
 *   n          = number of inner grid points per side (default 62, giving 64×64 grid)
 *   iterations = number of time steps              (default 10000)
 *
 * Domain:   unit square [0,1]×[0,1]
 * PDE:      u_t = c * (u_xx + u_yy)
 * BC:       left=10, right=40, front=30, back=50   (constant, as in paper §3)
 * IC:       all inner points = 0                   (as in paper §3)
 * Stability: dt = (ds^2) / (4*c)                  (as in paper §3)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

 /* ── boundary conditions (paper §3 example) ── */
static double bc_left = 10.0;   /* u(t, 0,   y) */
static double bc_right = 40.0;   /* u(t, 1,   y) */
static double bc_front = 30.0;   /* u(t, x,   0) */
static double bc_back = 50.0;   /* u(t, x,   1) */

/* ── initial condition for inner points ── */
static inline double f_init(int /*i*/, int /*j*/) { return 0.0; }

/* ── helpers ── */
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

    /* ── parse command-line args ── */
    int n = (argc > 1) ? atoi(argv[1]) : 62;   /* inner pts per side */
    int iterations = (argc > 2) ? atoi(argv[2]) : 10000;
    if (n < 1) n = 62;
    if (iterations < 1) iterations = 10000;

    /*
     * Grid layout (paper §2):
     *   Total points per side = N = n + 2  (0 … n+1, boundary inclusive)
     *   Spatial step ds = 1 / (n+1)
     *   Time step dt = ds^2 / (4*c)  — stability condition from paper §3
     */
    const int    N = n + 2;          /* total points per side */
    const double c = 0.1;            /* diffusivity (paper §3) */
    const double ds = 1.0 / (n + 1);  /* = 1/(N-1) */
    const double dt = (ds * ds) / (4.0 * c);
    const double r = c * dt / (ds * ds);   /* Fourier number; stable iff r <= 0.25 */

    /* ── domain decomposition along x (row) axis ── */
    /*
     * Each process owns a contiguous band of rows from the N×N grid.
     * Paper distributes along y; we distribute along i (rows), which is
     * equivalent by symmetry.
     *
     * Base share:  rows_per_proc = N / nprocs
     * Remainder:   extra = N % nprocs
     * Processes 0 … extra-1 each get one extra row.
     */
    int base_rows = N / nprocs;
    int extra = N % nprocs;

    /* number of owned rows and starting global row index for myid */
    int local_rows = base_rows + (myid < extra ? 1 : 0);
    int global_row0 = myid * base_rows + (myid < extra ? myid : extra);

    /*
     * Halo layout (paper §2):
     *   Local matrix has local_rows + 2 rows (one ghost row above and below).
     *   Ghost row 0        mirrors last row of process myid-1.
     *   Ghost row lr+1     mirrors first row of process myid+1.
     *   Owned rows: 1 … local_rows  (indices in local storage).
     */
    int lr = local_rows;          /* shorthand */
    int cols = N;                   /* all processes own every column */
    int local_total_rows = lr + 2;  /* including two ghost rows */

    double* u_old = alloc2d(local_total_rows, cols);
    double* u_new = alloc2d(local_total_rows, cols);

    /* ── initialise: boundary + initial conditions ── */
    /*
     * Global row g, column j:
     *   g == 0         → back   boundary   (bc_back)
     *   g == N-1       → front  boundary   (bc_front)
     *   j == 0         → left   boundary   (bc_left)
     *   j == N-1       → right  boundary   (bc_right)
     *   otherwise      → f_init(g, j) = 0
     */
    for (int li = 0; li <= lr + 1; li++) {
        int gi = global_row0 - 1 + li;   /* global row index for local row li */
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

    /*
     * ifirst / ilast (paper §2):
     *   First and last row index (in local storage) to be updated.
     *   Default: 1 … lr.
     *   First process: ifirst++ (its row 1 is the global top boundary, fixed).
     *   Last  process: ilast--  (its row lr is the global bottom boundary, fixed).
     */
    int ifirst = 1;
    int ilast = lr;
    if (myid == 0)          ifirst++;
    if (myid == nprocs - 1) ilast--;

    /* ── MPI derived datatype for a single row ── */
    MPI_Datatype row_type;
    MPI_Type_contiguous(cols, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    /* ── timing ── */
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    /* ── main time-stepping loop (paper §2) ── */
    for (int iter = 0; iter < iterations; iter++) {

        /* ── 1. compute interior points ── */
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

        /* ── 2. halo exchange (paper §2) ── */
        /*
         * Send my last owned row (u_new[lr]) to process myid+1's ghost row [0].
         * Receive into ghost row [lr+1] from process myid+1.
         * Send my first owned row (u_new[1]) to process myid-1's ghost row [lr+1].
         * Receive into ghost row [0] from process myid-1.
         *
         * Use non-blocking sends with blocking receives to avoid deadlock.
         */
        MPI_Request reqs[4];
        int nreqs = 0;

        if (myid < nprocs - 1) {
            /* send bottom owned row to right neighbour */
            MPI_Isend(&u_new[IDX(lr, 0, cols)], 1, row_type,
                myid + 1, 0, MPI_COMM_WORLD, &reqs[nreqs++]);
            /* receive into bottom ghost row from right neighbour */
            MPI_Irecv(&u_new[IDX(lr + 1, 0, cols)], 1, row_type,
                myid + 1, 1, MPI_COMM_WORLD, &reqs[nreqs++]);
        }
        if (myid > 0) {
            /* send top owned row to left neighbour */
            MPI_Isend(&u_new[IDX(1, 0, cols)], 1, row_type,
                myid - 1, 1, MPI_COMM_WORLD, &reqs[nreqs++]);
            /* receive into top ghost row from left neighbour */
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

    /* ── gather result on rank 0 and print summary ── */
    /*
     * Each process sends its local_rows owned rows (li = 1…lr) to rank 0.
     * Rank 0 assembles the full N×N grid and reports a few diagnostics.
     */

     /* build per-process send counts and displacements (in rows) */
    int* rcounts = NULL, * displs = NULL;
    double* full_grid = NULL;

    if (myid == 0) {
        rcounts = (int*)malloc(nprocs * sizeof(int));
        displs = (int*)malloc(nprocs * sizeof(int));
        full_grid = alloc2d(N, cols);   /* N × N */

        int offset = 0;
        for (int p = 0; p < nprocs; p++) {
            int pr = base_rows + (p < extra ? 1 : 0);
            /* skip boundary rows for rank 0 and last rank */
            int owned;
            if (p == 0 && nprocs == 1) owned = pr - 2;  /* both boundaries */
            else if (p == 0 || p == nprocs - 1) owned = pr - 1;
            else                                owned = pr;
            rcounts[p] = owned * cols;
            displs[p] = offset * cols;
            offset += owned;
        }
    }

    /* pack local owned rows into a contiguous send buffer */
    int owned_rows;
    if (nprocs == 1)          owned_rows = lr - 2;  /* both bdy */
    else if (myid == 0 || myid == nprocs - 1) owned_rows = lr - 1;
    else                           owned_rows = lr;

    int send_start = ifirst;   /* first row to send (== ifirst for all) */
    double* send_buf = &u_old[IDX(send_start, 0, cols)];

    MPI_Gatherv(send_buf, owned_rows * cols, MPI_DOUBLE,
        full_grid, rcounts, displs, MPI_DOUBLE,
        0, MPI_COMM_WORLD);

    /* ── rank 0: print results ── */
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

        /* Sample temperatures at a few interior points */
        printf("  Sample temperatures (inner grid):\n");
        int pts[][2] = { {N / 4, N / 4}, {N / 2, N / 2}, {3 * N / 4, 3 * N / 4},
                        {N / 4, 3 * N / 4}, {3 * N / 4, N / 4} };
        int npts = (int)(sizeof(pts) / sizeof(pts[0]));
        for (int k = 0; k < npts; k++) {
            int gi = pts[k][0], gj = pts[k][1];
            /* full_grid stores only inner rows (ifirst_global…ilast_global) */
            /* adjust index: rank 0 skipped row 0, so full_grid row 0 = global row 1 */
            int fi = gi - 1;
            if (fi >= 0 && fi < N - 2)
                printf("    u(x=%d/%d, y=%d/%d) = %.4f\n",
                    gi, N - 1, gj, N - 1,
                    full_grid[IDX(fi, gj, cols)]);
        }

        /* min / max over interior */
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

    /* ── report per-process timing via reduce ── */
    double max_time;
    MPI_Reduce(&elapsed, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myid == 0)
        printf("  Max wall time across all ranks: %.6f s\n\n", max_time);

    /* ── cleanup ── */
    MPI_Type_free(&row_type);
    free(u_old); free(u_new);

    MPI_Finalize();
    return 0;
}