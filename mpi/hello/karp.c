#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#define f(x) ((double)(4.0/(1.0+x*x)))
#define pi ((double)(4.0*atan(1.0)))

main(int argc, char **argv)

{

  /* This simple program approximates pi by computing pi = integral
   * from 0 to 1 of 4/(1+x*x)dx which is approximated by sum 
   * from k=1 to N of 4 / (1+((k-.5)/N)**2).  The only input data 
   * required is N.                                       
   * 
   * 3/21/97 RLF Change floats to doubles
   * 8/97 SHM Read input from file to accommodate VW Companion
   */

  /* Each process is given a chunk of the interval to do. */

  double err, sum, w, tmpsum;
  int i, N, int_start, int_end;
  FILE *GetInterval;

  int rank, size, tag, rc;
  MPI_Status status;


  /* Insert call to startup routine that returns the number of tasks and
   * the taskid of the current instance.
   */

  rc = MPI_Init(&argc, &argv);
  rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
  rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  tag = 100;

  /* 
   * Now solicit a new value for N.  When it is 0, then you should depart.
   * This would be a good place to unenroll yourself as well.
   */



  if (rank == 0) { //Master
    // Reads the file and distributes the pieces among workers
    GetInterval = fopen( "./data_values", "r" );
    fscanf(GetInterval, "%d", &N);
    printf("Approximation interval is %d\n", N);
   }
   rc = MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

  while (N > 0) 
  {
    int_start = 1 + (int)(N/size) * rank;
    int_end = (int)(N/size) * (rank+1);
    if (rank == (size-1))
      int_end = N;

    w = 1.0/(double)N;
    sum = 0.0;
    //for (i = 1; i <= N; i++)
    for (i = int_start; i <= int_end; i++)
      sum = sum + f(((double)i-0.5)*w);
    sum = sum * w;

    MPI_Reduce(&sum, &tmpsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      err = tmpsum - pi;
      printf("sum, err = %7.5f, %10e\n", tmpsum, err);
      fscanf(GetInterval, "%d", &N);
      printf("Approximation interval is %d\n", N);
      rc = MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
      rc = MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }
  if (rank == 0)
    fclose( GetInterval );
  MPI_Finalize();
}
