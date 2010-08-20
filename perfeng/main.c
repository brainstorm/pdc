#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define BLOCK 30 
#define min(a, b) ((a < b) ? a : b)
#define MUL_OUTPUT 0
//TODO: sse

//#ifdef GNU
	#include <xmmintrin.h>
//#else
//Intel Compiler
	#include <emmintrin.h>
	#include <mmintrin.h>
//#endif

// example prototype for your matmul function
void mul(double* dest, const double* a, const double* b, int N);

inline double gettime(void) 
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec + 1e-6*tv.tv_usec;
}

inline void fill_mat(double* m, int M, double val)
{
    int i, j;
    for (i=0; i<M; i++)
	for (j=0; j<M; j++)
	    m[M*i + j] = val;
}

inline void mul(double* dest, const double* a, const double* b, int M){
    int i, j, k, jj, kk;
    for (jj=0; jj<M; jj+=BLOCK) 
	for (kk=0; kk<M; kk+=BLOCK)
	    for (i=0; i<M; ++i)
		for (j=jj; j<min(jj+BLOCK, M); ++j){
		    double sum = 0.0;
		    for (k=kk; k<min(kk+BLOCK, M); ++k)
			sum += a[M*i + k] * b[M*k + j];
		    dest[M*i + j] += sum;
		}
}


void mul_sse(double* dest, const double* a, const double* b, int M){
	int i, j, k, l;

	__m128d *ae, *be, *res;

//	for (i=0; i<M; i++) {
		ae = (__m128d*) &a[0];
		be = (__m128d*) &b[0];
		res = (__m128d*) &dest[0];

		*res = __mm_mul_pd(*ae, *be);
//__mm_malloc
//	}

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
    if (MUL_OUTPUT) {    
	for (i=0; i<M; i++, printf("\n"))
	    for (j=0; j<M; j++, printf(" "))
		printf("%lf", C[M*i + j]);
    }
    free(A);
    free(B);
    free(C);
    return 0;
}
