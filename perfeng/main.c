#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

// example prototype for your matmul function
void mul(double* dest, const double* a, const double* b, int N);

double gettime(void) 
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec + 1e-6*tv.tv_usec;
}

void fill_mat(double* m, int M, double val)
{
    for(int i=0; i<M; i++)
	for(int j=0; j<M; j++)
	    m[M*i + j] = val;
}

int main(int args, char* argv[])
{
    const int M = atoi(argv[1]);

    double* A = (double*) malloc(M*M*sizeof(double));
    double* B = (double*) malloc(M*M*sizeof(double));
    double* C = (double*) malloc(M*M*sizeof(double));

    fill_mat(A,M,1.0);
    fill_mat(B,M,2.0);
    fill_mat(C,M,0.0);

    double t = gettime();
    mul(C,A,B,M);
    t = gettime()-t;

    printf("%d\t%f\t%E\n",M,t,2*pow(M,3)/t);

    free(A);
    free(B);
    free(C);
    return 0;
}
