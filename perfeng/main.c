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
#define IFSIZE 32

#define min(a, b) ((a < b) ? a : b)
#define max(a, b) ((a > b) ? a : b)
#define MUL_OUTPUT 0
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
  __m128d ae1, be1, ae2, be2, ae3, be3;
  __m128d sum1, sum2, sum3;
  __m128d res1, res2, res3;
  if (M>(IFSIZE*BLOCK2)) {
#pragma omp parallel for default(none) shared(a, b, dest, M) private(i0, j0, k0, i1, i, j1, j, k1, k, sum, ae, be, res, dummy, res1, res2, res3, be1, be2, be3, sum1, sum2, sum3) schedule(dynamic)
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
                  for (j=j1; j<min(j1+BLOCK1, M); j+=8){ 
                    sum = _mm_setzero_pd();
                    sum1 = _mm_setzero_pd();
                    sum2 = _mm_setzero_pd();
                    sum3 = _mm_setzero_pd();
                    for (k=k1; k<min(k1+BLOCK1, M); ++k) {

                      // Loading values into __m128d
                      ae = _mm_load1_pd(&(a[M*i+k]));

                      be = _mm_load_pd(&(b[M*k+j])); 
                      be1 = _mm_load_pd(&(b[M*k+j+2])); 
                      be2 = _mm_load_pd(&(b[M*k+j+4])); 
                      be3 = _mm_load_pd(&(b[M*k+j+6])); 


                      // Performing multiplication and add (sum += a * b)
                      sum = _mm_add_pd(sum, _mm_mul_pd(ae, be)); 
                      sum1 = _mm_add_pd(sum1, _mm_mul_pd(ae, be1)); 
                      sum2 = _mm_add_pd(sum2, _mm_mul_pd(ae, be2)); 
                      sum3 = _mm_add_pd(sum3, _mm_mul_pd(ae, be3)); 

                      //_mm_store_pd(dummy, sum);
                      //printf("dummy: %lf %lf\n", dummy[0], dummy[0]);
                    }
                    //printf("dest: %lf %lf\n", dest[M*i+j], dest[M*i+j+1]);

                    // Add result
                    res = _mm_load_pd(&(dest[M*i+j]));
                    res1 = _mm_load_pd(&(dest[M*i+j+2]));
                    res2 = _mm_load_pd(&(dest[M*i+j+4]));
                    res3 = _mm_load_pd(&(dest[M*i+j+6]));

                    res = _mm_add_pd(res, sum);
                    res1 = _mm_add_pd(res1, sum1);
                    res2 = _mm_add_pd(res2, sum2);
                    res3 = _mm_add_pd(res3, sum3);

                    _mm_store_pd(&(dest[M*i+j]), res);
                    _mm_store_pd(&(dest[M*i+j+2]), res1);
                    _mm_store_pd(&(dest[M*i+j+4]), res2);
                    _mm_store_pd(&(dest[M*i+j+6]), res3);

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
#pragma omp parallel for default(none) shared(a, b, dest, M) private(i1, i, j1, j, k1, k, sum, ae, be, res, dummy, res1, res2, res3, be1, be2, be3, sum1, sum2, sum3) schedule(dynamic)
    // L1 block
    for (i1=0; i1<M; i1+=BLOCK1) {
      for (j1=0; j1<M; j1+=BLOCK1) {
        for (k1=0; k1<M; k1+=BLOCK1) {
          // Multiplication loop
          for (i=i1; i<min(i1+BLOCK1, M); ++i) {
            for (j=j1; j<min(j1+BLOCK1, M); j+=8){ 
              sum = _mm_setzero_pd();
              sum1 = _mm_setzero_pd();
              sum2 = _mm_setzero_pd();
              sum3 = _mm_setzero_pd();
              for (k=k1; k<min(k1+BLOCK1, M); ++k) {
                // Loading values into __m128d
                ae = _mm_load1_pd(&(a[M*i+k]));
                
                be = _mm_load_pd(&(b[M*k+j])); 
                be1 = _mm_load_pd(&(b[M*k+j+2])); 
                be2 = _mm_load_pd(&(b[M*k+j+4])); 
                be3 = _mm_load_pd(&(b[M*k+j+6])); 

                // Performing multiplication and add (sum += a * b)
                sum = _mm_add_pd(sum, _mm_mul_pd(ae, be)); 
                sum1 = _mm_add_pd(sum1, _mm_mul_pd(ae, be1)); 
                sum2 = _mm_add_pd(sum2, _mm_mul_pd(ae, be2)); 
                sum3 = _mm_add_pd(sum3, _mm_mul_pd(ae, be3)); 

                // Performing multiplication and add (sum += a * b)
                //sum = _mm_add_pd(sum, _mm_mul_pd(ae, be)); 
                //_mm_store_pd(dummy, sum);
              }
              // Add result
              res = _mm_load_pd(&(dest[M*i+j]));
              res1 = _mm_load_pd(&(dest[M*i+j+2]));
              res2 = _mm_load_pd(&(dest[M*i+j+4]));
              res3 = _mm_load_pd(&(dest[M*i+j+6]));

              res = _mm_add_pd(res, sum);
              res1 = _mm_add_pd(res1, sum1);
              res2 = _mm_add_pd(res2, sum2);
              res3 = _mm_add_pd(res3, sum3);

              _mm_store_pd(&(dest[M*i+j]), res);
              _mm_store_pd(&(dest[M*i+j+2]), res1);
              _mm_store_pd(&(dest[M*i+j+4]), res2);
              _mm_store_pd(&(dest[M*i+j+6]), res3);
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
  FILE *sse_file;
  sse_file = fopen("/scratch/sse_corr.out", "w");
  if (sse_file != NULL) {
    for (i=0; i<M; i++, fprintf(sse_file, "\n"))
      for (j=0; j<M; j++, fprintf(sse_file, " "))
        fprintf(sse_file, "%lf", C[M*i + j]);
    fclose(sse_file);
  }
#endif        
  
#ifdef USE_BLAS
  fill_mat(C,M,0.0);
  t = gettime();
  mul_blas(C,A,B,M);
  t = gettime()-t;

  printf("%d\t%f\t%E\n",M,t,2*pow(M,3)/t);
  
#ifdef TEST_CORRECTNESS
  FILE *blas_file;
  blas_file = fopen("/scratch/blas_corr.out", "w");
  if (blas_file != NULL) {
    for (i=0; i<M; i++, fprintf(blas_file, "\n"))
      for (j=0; j<M; j++, fprintf(blas_file, " "))
        fprintf(blas_file, "%lf", C[M*i + j]);
    fclose(blas_file);
  }
#endif
#endif  
  
  _mm_free(A);
  _mm_free(B);
  _mm_free(C);
  return 0;
}
