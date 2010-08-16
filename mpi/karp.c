#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define f(x) ((double)(4.0/(1.0+x*x)))
#define pi ((double)(4.0*atan(1.0)))

main()

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

    double err, sum, w;
    int i, N;
    FILE *GetInterval;


    /* Insert call to startup routine that returns the number of tasks and
     * the taskid of the current instance.
     */

    /* 
     * Now solicit a new value for N.  When it is 0, then you should depart.
     * This would be a good place to unenroll yourself as well.
     */

    GetInterval = fopen( "./data_values", "r" );
    fscanf(GetInterval, "%d", &N);
    printf("Approximation interval is %d\n", N);

    while (N > 0) 
    {
	w = 1.0/(double)N;
	sum = 0.0;
	for (i = 1; i <= N; i++)
	    sum = sum + f(((double)i-0.5)*w);
	sum = sum * w;
	err = sum - pi;
	printf("sum, err = %7.5f, %10e\n", sum, err);
	fscanf(GetInterval, "%d", &N);
	printf("Approximation interval is %d\n", N);
    }
    fclose( GetInterval );
}
