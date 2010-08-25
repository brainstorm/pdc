#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>

#ifdef USE_BLAS
#include "mkl_cblas.h"
#endif

// Bruteforcing...
#define BLOCK2 256
#define BLOCK1 32
#define IFSIZE 2

#define min(a, b) ((a < b) ? a : b)
#define max(a, b) ((a > b) ? a : b)
#define MUL_OUTPUT 1
//TODO: sse, profile value of BLOCK2 (bruteforce) or maybe use cache line size somehow
// omp inner loop ? profile where time is spent

/*
According to /proc/cpuinfo:
(...)
cache_alignment : 64

Meaning: 64 *bytes* cache line

It would be reasonable then to set BLOCK to 64 though

*/

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
      m[M*i + j] = val*(rand()/(double)(RAND_MAX+1.0));
      //m[M*i + j] = val*(rand()%10);
}

#ifdef USE_BLAS
void mul_blas(double* c, double* a, double* b, int N){
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
#endif


inline void mul(double* dest, const double* a, const double* b, int M){
  int i, j, k, jj, kk;

  for (jj=0; jj<M; jj+=BLOCK2) 
    for (kk=0; kk<M; kk+=BLOCK2)
      for (i=0; i<M; ++i)
        for (j=jj; j<min(jj+BLOCK2, M); ++j){
          double sum = 0.0;
          for (k=kk; k<min(kk+BLOCK2, M); ++k)
            sum += a[M*i + k] * b[M*k + j];
          dest[M*i + j] += sum;
        }
}

inline void mul_trans(double* dest, const double* a, const double* b, int M){
  int i, j, k, jj, kk;
  
  // Transposing matrix b
  double* bT = (double*) malloc(M*M*sizeof(double));
  for (i=0; i<M; ++i)
    for (j=0; j<M; ++j)
      bT[M*i+j] = b[M*j+i];


  for (jj=0; jj<M; jj+=BLOCK2) 
    for (kk=0; kk<M; kk+=BLOCK2)
      for (i=0; i<M; ++i)
        for (j=jj; j<min(jj+BLOCK2, M); ++j){
          double sum = 0.0;
          for (k=kk; k<min(kk+BLOCK2, M); ++k)
            sum += a[M*i + k] * bT[M*j + k];
          dest[M*i + j] += sum;
        }

   //XXX
  //_mm_free(bT);      
}


void mul_sse(double* restrict dest, const double* restrict a, const double* restrict b, int M){
  int i0, i1, i, j0, j1, j, k0, k1, k;
  double dummy[2];

  __m128d ae, be, res, sum, tmp;
  if (M>(IFSIZE*BLOCK2)) {
#pragma omp parallel for default(none) shared(a, b, dest, M) private(i0, i1, i, j0, j1, j, k0, k1, k, sum, ae, be, res, dummy) schedule(dynamic)
    // L2 block
    for (i0=0; i0<M; i0+=BLOCK2) {
      for (j0=0; j0<M; j0+=BLOCK2) {
        for (k0=0; k0<M; k0+=BLOCK2) {
          // L1 block
          for (i1=i0; i1<min(i0+BLOCK2, M); i1+=BLOCK1) {
            for (j1=j0; j1<min(j0+BLOCK2, M); j1+=BLOCK1) {
              for (k1=k0; k1<min(k0+BLOCK2, M); k1+=BLOCK1) {
                // Multiplication loop
                for (i=i1; i<min(i1+BLOCK1, M); ++i) {
                  for (j=j1; j<min(j1+BLOCK1, M); j+=2){ 
                    sum = _mm_setzero_pd();
                    for (k=k1; k<min(k1+BLOCK1, M); ++k) {

                      // Loading values into __m128d
                      ae = _mm_load1_pd(&(a[M*i+k]));
                      be = _mm_load_pd(&(b[M*k+j])); 

                      // Performing multiplication and add (sum += a * b)
                      sum = _mm_add_pd(sum, _mm_mul_pd(ae, be)); 

                      //_mm_store_pd(dummy, sum);
                      //printf("dummy: %lf %lf\n", dummy[0], dummy[0]);
                    }
                    //printf("dest: %lf %lf\n", dest[M*i+j], dest[M*i+j+1]);
                    // Add result
                    res = _mm_load_pd(&(dest[M*i+j]));
                    res = _mm_add_pd(res, sum);
                    _mm_store_pd(&(dest[M*i+j]), res);

                    //_mm_store_pd(dummy, sum);
                    //printf("dummy: %lf %lf\n", dummy[0], dummy[1]);
                    //printf("dest: %lf %lf\n", dest[0], dest[1]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  else {
    // L1 block
    for (i1=0; i1<M; i1+=BLOCK1) {
      for (j1=0; j1<M; j1+=BLOCK1) {
        for (k1=0; k1<M; k1+=BLOCK1) {
#pragma omp parallel for default(none) shared(a, b, dest, M) private(i1, i, j1, j, k1, k, sum, ae, be, res, dummy) schedule(dynamic)
          // Multiplication loop
          for (i=i1; i<min(i1+BLOCK1, M); ++i) {
            for (j=j1; j<min(j1+BLOCK1, M); j+=2){ 
              sum = _mm_setzero_pd();
              for (k=k1; k<min(k1+BLOCK1, M); ++k) {
                // Loading values into __m128d
                ae = _mm_load1_pd(&(a[M*i+k]));
                be = _mm_load_pd(&(b[M*k+j])); 
                
                // Performing multiplication and add (sum += a * b)
                sum = _mm_add_pd(sum, _mm_mul_pd(ae, be)); 
                _mm_store_pd(dummy, sum);
              }
              // Add result
              res = _mm_load_pd(&(dest[M*i+j]));
              res = _mm_add_pd(res, sum);
              _mm_store_pd(&(dest[M*i+j]), res);
            }
          }
        }
      }
    }
  }

}

int main(int args, char* argv[])
{
  int i, j;
  const int M = atoi(argv[1]);
  srand(5);

  double* A = (double*) _mm_malloc(M*M*sizeof(double), 16);
  double* B = (double*) _mm_malloc(M*M*sizeof(double), 16);
  double* C = (double*) _mm_malloc(M*M*sizeof(double), 16);

  fill_mat(A,M,1.0);
  fill_mat(B,M,2.0);
  fill_mat(C,M,0.0);
  
#ifdef TEST_CORRECTNESS
  printf("Matrix A\n");
    for (i=0; i<M; i++, printf("\n"))
      for (j=0; j<M; j++, printf(" "))
        printf("%lf", A[M*i + j]);
  
  printf("Matrix B\n");
    for (i=0; i<M; i++, printf("\n"))
      for (j=0; j<M; j++, printf(" "))
        printf("%lf", B[M*i + j]);
#endif
  
  double t = gettime();
  //mul(C,A,B,M);
  mul_sse(C,A,B,M);
  t = gettime()-t;

  printf("%d\t%f\t%E\n",M,t,2*pow(M,3)/t);
  
#ifdef TEST_CORRECTNESS
  printf("SSE output\n");  
    for (i=0; i<M; i++, printf("\n"))
      for (j=0; j<M; j++, printf(" "))
        printf("%lf", C[M*i + j]);
#endif        
  
#ifdef USE_BLAS
  fill_mat(C,M,0.0);
  t = gettime();
  mul_blas(C,A,B,M);
  t = gettime()-t;

  printf("%d\t%f\t%E\n",M,t,2*pow(M,3)/t);
  
#ifdef TEST_CORRECTNESS
  printf("BLAS output\n");
    for (i=0; i<M; i++, printf("\n"))
      for (j=0; j<M; j++, printf(" "))
        printf("%lf", C[M*i + j]);
#endif
#endif  
  
  _mm_free(A);
  _mm_free(B);
  _mm_free(C);
  return 0;
}
