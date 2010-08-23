#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>

#define BLOCK 512
#define BLOCK2 128
#define min(a, b) ((a < b) ? a : b)
#define MUL_OUTPUT 0
//TODO: sse, profile value of BLOCK (bruteforce) or maybe use cache line size somehow

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

inline void mul_trans(double* dest, const double* a, const double* b, int M){
  int i, j, k, jj, kk;
  
  // Transposing matrix b
  double* bT = (double*) malloc(M*M*sizeof(double));
  for (i=0; i<M; ++i)
    for (j=0; j<M; ++j)
      bT[M*i+j] = b[M*j+i];


  for (jj=0; jj<M; jj+=BLOCK) 
    for (kk=0; kk<M; kk+=BLOCK)
      for (i=0; i<M; ++i)
        for (j=jj; j<min(jj+BLOCK, M); ++j){
          double sum = 0.0;
          for (k=kk; k<min(kk+BLOCK, M); ++k)
            sum += a[M*i + k] * bT[M*j + k];
          dest[M*i + j] += sum;
        }

  _mm_free(bT);      
}


void mul_sse(double* restrict dest, const double* restrict a, const double* restrict b, int M){
  int i0, i1, i, j0, j1, j, k0, k1, k;
  double dummy[2];
  // Transposing matrix b
  /*double* bT = (double*) _mm_malloc(M*M*sizeof(double), 16);
  for (i=0; i<M; ++i)
    for (j=0; j<M; ++j)
      bT[M*i+j] = b[M*j+i];*/

  __m128d ae, be, res, sum;

#pragma omp parallel for default(none) shared(a, b, dest) private(i0, i1, i, j0, j1, j, k0, k1, k, M, sum, ae, be, res, dummy)
  // First block
  for (i0=0; i0<M; i0+=BLOCK) {
    for (j0=0; j0<M; j0+=BLOCK) {
      for (k0=0; k0<M; k0+=BLOCK) {
        // Second block
        for (i1=i0; i1<i0+BLOCK; i1+=BLOCK2) {
          for (j1=j0; j1<j0+BLOCK; j1+=BLOCK2) {
            for (k1=k0; k1<k0+BLOCK; k1+=BLOCK2) {
              // Multiplication loop
              for (i=i1; i<min(i1+BLOCK2, M); ++i) {
                res = _mm_setzero_pd();
                for (j=j1; j<min(j1+BLOCK2, M); j+=2){ 
                  for (k=k1; k<min(k1+BLOCK2, M); k+=2) {
                    sum = _mm_setzero_pd();

                    // Loading values into __m128d
                    ae = _mm_load1_pd(&(a[M*i+k]));
                    be = _mm_load_pd(&(b[M*k+j])); 

                    // Performing multiplication and add (sum += a * b)
                    sum = _mm_add_pd(sum, _mm_mul_pd(ae, be)); 

                    //_mm_store_pd(dummy, sum);
                    //printf("dummy: %lf %lf\n", dummy[0], dummy[0]);
                  }
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

  //_mm_free(bT);    
}

int main(int args, char* argv[])
{
  int i, j;
  const int M = atoi(argv[1]);

  double* A = (double*) _mm_malloc(M*M*sizeof(double), 16);
  double* B = (double*) _mm_malloc(M*M*sizeof(double), 16);
  double* C = (double*) _mm_malloc(M*M*sizeof(double), 16);

  fill_mat(A,M,1.0);
  fill_mat(B,M,2.0);
  fill_mat(C,M,0.0);
  
  double t = gettime();
  //mul(C,A,B,M);
  mul_sse(C,A,B,M);
  t = gettime()-t;

  printf("%d\t%f\t%E\n",M,t,2*pow(M,3)/t);
  
  fill_mat(C,M,0.0);
  t = gettime();
  mul_trans(C,A,B,M);
  //mul_sse(C,A,B,M);
  t = gettime()-t;

  printf("%d\t%f\t%E\n",M,t,2*pow(M,3)/t);
  
  
  if (MUL_OUTPUT) {    
    for (i=0; i<M; i++, printf("\n"))
      for (j=0; j<M; j++, printf(" "))
        printf("%lf", C[M*i + j]);
  }
  _mm_free(A);
  _mm_free(B);
  _mm_free(C);
  return 0;
}
