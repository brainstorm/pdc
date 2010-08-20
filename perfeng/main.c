#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#ifdef GNU
	#include <xmmintrin.h>
#else
//Intel Compiler
	#include <emmintrin.h>
#endif

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
    int i, j;
    for (i=0; i<M; i++)
	for (j=0; j<M; j++)
	    m[M*i + j] = val;
}

void mul(double* dest, const double* a, const double* b, int M){
	int i, j, k, l;
	for (i=0; i<M; i++)
		for (j=0; j<M; j++){
			double sum = 0.0;
			for (k=0; k<M; k++)
				sum += a[M*i + k] * b[M*k + j];
			dest[M*i + j] = sum;
		}
}


void mul_sse(double* dest, const double* a, const double* b, int M){
	int i, j, k, l;

	__m128d *ae, *be, *res;

	for (i=0; i<M; i++) {
		ae = (_m128d*) &a[i];
		be = (_m128d*) &b[i];
		*res = _mm_mul_sd(ae, be);
	}

//_mm_madd_epi16

/*
	for (i=0; i<M; i++)
		for (j=0; j<M; j++){
			double sum = 0.0;
			for (k=0; k<M; k++)
				sum += a[M*i + k] * b[M*k + j];
			dest[M*i + j] = sum;
		}
*/
}

int main(int args, char* argv[])
{
    int i, j;
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
/*
    for (i=0; i<M; i++, printf("\n"))
	for (j=0; j<M; j++, printf(" "))
		printf("%lf", C[M*i + j]);
*/
    free(A);
    free(B);
    free(C);
    return 0;
}
