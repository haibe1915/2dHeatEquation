//#define MSMPI_NO_SAL
//#define _Out_writes_(x)
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <mpi.h>
//
//#define  M       20
//#define  dx      0.1
//#define  tol     0.000001
//
//void KhoiTao(float* Told) {
//    int i;
//    for (i = 0; i < M; i++)
//        *(Told + i) = 25.0;
//}
//
//
//
//int main(int argc, char* argv[]) {
//    int i;
//    float T_left = 100.0, T_right = 25.0;
//    float* T_old, * T_new;
//
//    T_old = (float*)malloc((M + 1) * sizeof(float));
//    T_new = (float*)malloc((M + 1) * sizeof(float));
//
//    KhoiTao(T_old);
//
//    int Rank, NP;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &NP);
//
//    int Mc = M / NP;
//
//    float* T_old_c = (float*)malloc((Mc + 2) * sizeof(float));
//    float* T_new_c = (float*)malloc((Mc + 2) * sizeof(float));
//
//    for (int i = 1; i <= Mc; i++) {
//        T_old_c[i] = T_right;
//    }
//
//    if (Rank == 0) T_old_c[1] = T_left;
//    if (Rank == NP - 1) T_new_c[Mc] = T_right;
//
//    int iter = 0;
//    float max_diff, global_diff;
//
//    do {
//        iter++;
//        max_diff = 0.0;
//
//        MPI_Status status;
//
//        if (Rank < NP - 1)
//            MPI_Send(&T_old_c[Mc], 1, MPI_FLOAT, Rank + 1, 0, MPI_COMM_WORLD);
//        if (Rank > 0)
//            MPI_Recv(&T_old_c[0], 1, MPI_FLOAT, Rank - 1, 0, MPI_COMM_WORLD, &status);
//
//        if (Rank > 0)
//            MPI_Send(&T_old_c[1], 1, MPI_FLOAT, Rank - 1, 1, MPI_COMM_WORLD);
//        if (Rank < NP - 1)
//            MPI_Recv(&T_old_c[Mc + 1], 1, MPI_FLOAT, Rank + 1, 1, MPI_COMM_WORLD, &status);
//
//        if (Rank == 0) T_old_c[0] = T_left;
//        if (Rank == NP - 1) T_old_c[Mc + 1] = T_right;
//
//
//        for (int i = 1; i <= Mc; i++) {
//            T_new_c[i] = (T_old_c[i - 1] + T_old_c[i + 1]) / 2.0;
//
//            float diff = fabs(T_new_c[i] - T_old_c[i]);
//            if (diff > max_diff) max_diff = diff;
//        }
//
//        for (int i = 1; i <= Mc; i++) {
//            T_old_c[i] = T_new_c[i];
//        }
//
//        MPI_Allreduce(&max_diff, &global_diff, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
//
//    } while (global_diff > tol);
//
//
//    MPI_Gather(&T_old_c[1], Mc, MPI_FLOAT,
//        T_old, Mc, MPI_FLOAT,
//        0, MPI_COMM_WORLD);
//
//    if (Rank == 0) {
//        printf("Iterations: %d\n", iter);
//        T_old[0] = T_left;
//        T_old[M] = T_right;
//
//        for (int i = 0; i <= M; i++) {
//            printf("%.6f ", T_old[i]);
//            if (i % 10 == 0) printf("\n");
//        }
//    }
//
//    MPI_Finalize();
//    return 0;
//}