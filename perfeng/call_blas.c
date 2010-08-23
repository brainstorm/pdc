#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include "mkl_cblas.h"

#define OUTPUT 0
/*
SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      ..
      .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*/

void mul(double* c, double* a, double* b, int N){
 /* parameters for C <- A*B, N-by-N */
    const CBLAS_ORDER order = CblasRowMajor;
    const CBLAS_TRANSPOSE transA = CblasNoTrans;
    const CBLAS_TRANSPOSE transB = CblasNoTrans;
    const double alpha = 1.0;
    const double beta = 0.0;
    const int lda = N;
    const int ldb = N;
    const int ldc = N;

    cblas_dgemm(order, transA, transB, N, N, N, alpha,
		a, lda, b, ldb, beta, c, ldc);
}

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
  if (OUTPUT) {
    for (i=0; i<M; i++, printf("\n"))
      for (j=0; j<M; j++, printf(" "))
        printf("%lf", C[M*i + j]);
  }	

  free(A);
  free(B);
  free(C);
  return 0;
}
