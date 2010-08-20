#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

double gettime(void) 
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec + 1e-6*tv.tv_usec;
}

int main(int args, char* argv[])
{
    const int N = atoi(argv[1]);
    double* mat = (double*) malloc(N*N*sizeof(double));
    int i,j;

    for(i=0; i<N*N; i++) /* initialize some values */
        mat[i] = i;

    double t = gettime();

     for(j=0; j<N; j++)
   for(i=0; i<N; i++) 
            mat[i*N + j] *= 10.2; /* scale each value*/

    t = gettime()-t;
    printf("N = %d\t time: %f seconds \t flops: %e\n",N,t,(N/t)*N);

    free(mat);
}
