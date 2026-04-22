//#define MSMPI_NO_SAL
//#include <stdio.h>
//#define _Out_writes_(x)
//#include <mpi.h>
//
//int main(int argc, char** argv) {
//    int rank, size;
//    MPI_Status stat;
//
//    MPI_Init(&argc, &argv);
//
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    printf("Hello from process %d of %d\n", rank, size);
//    fflush(stdout);  /* <-- B?T BU?C khi ch?y qua mpiexec */
//
//    int s[10] = { 0,1,2,3,4,5,6,7,8,9 };
//    int r[10];
//    int dest = (rank + 1) % size;
//    int src = (rank - 1 + size) % size;
//
//    MPI_Sendrecv(
//        s, 10, MPI_INT, dest, 0,
//        r, 10, MPI_INT, src, 0,
//        MPI_COMM_WORLD, &stat
//    );
//
//    printf("Process %d received from %d: ", rank, src);
//    int i;
//    for (i = 0; i < 10; i++) {
//        printf("%d ", r[i]);
//    }
//    printf("\n");
//    fflush(stdout);  /* <-- flush l?n n?a d? ch?c ch?n */
//
//    MPI_Finalize();
//    return 0;
//}