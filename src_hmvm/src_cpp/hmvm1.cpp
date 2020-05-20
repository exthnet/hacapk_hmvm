#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hacapk.h"
#include "hmvm_seq.h"
#include "hmvm_omp.h"
#include "hmvm_cuda.h"
#include "hmvm_magma.h"

int main(int argc, char **argv)
{
  int i, ret;
  FILE *F;
  matrix<double> *matD=NULL;
  matrix2<double> *matD2=NULL;
  matrix<float> *matF=NULL;
  matrix2<float> *matF2=NULL;
  int dump_result=0;
  int nbench=5;

  if(argc<2){
	printf("usage: %s filename [dump(0,1) nbench]\n", argv[0]);
	return -1;
  }
  if(argc>=3){dump_result=atoi(argv[2]);printf("dump_result = %d\n",dump_result);}
  if(argc>=4){nbench=atoi(argv[3]);printf("nbench = %d\n",nbench);}
  printf("call loadmatrix %s\n", argv[1]); fflush(stdout);
  matD = new matrix<double>();
  matD2 = new matrix2<double>();
  ret = loadHmatrix(argv[1], matD, matD2);
  if(ret)return -1;

  int nd = matD->nd;
  printf("nd=%d\n",nd);
  double *b=NULL;
  b = new double[nd];//(double*)malloc(sizeof(double)*matD.nd);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
  }

  F = fopen("vec.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", b[i]);
  fclose(F);

  // test
  //for(i=0;i<96;i++)
  //i=0; hmvm_cuda1<double>(matD2, b, i, dump_result, nbench);

  // ######## ######## ######## ######## ######## ######## ######## ########
  // ready for hmvm
  // ######## ######## ######## ######## ######## ######## ######## ########

#if 1
  // ######## ######## ######## ######## ######## ######## ######## ########
  // double
  // ######## ######## ######## ######## ######## ######## ######## ########
  printf("######## ######## ######## ########\n");
  printf("double\n");
  printf("######## ######## ######## ########\n");
  // sequential
  hmvm_seq<double>(matD, matD2, b, dump_result, nbench);
  // OpenMP
  hmvm_omp<double>(matD, matD2, b, dump_result, nbench);
  // MKL
#ifdef __INTEL_COMPILER
  hmvm_blas_p(matD, matD2, b, dump_result, nbench);
  hmvm_blas_p_bench(matD, matD2, b, nbench);
  hmvm_blas_s(matD, matD2, b, dump_result, nbench);
  hmvm_blas_s_bench(matD, matD2, b, nbench);
  hmvm_cblas_p(matD, matD2, b, dump_result, nbench);
  hmvm_cblas_p_bench(matD, matD2, b, nbench);
  hmvm_cblas_s(matD, matD2, b, dump_result, nbench);
  hmvm_cblas_s_bench(matD, matD2, b, nbench);
#endif
  // CUDA
#ifdef _USE_CUDA
  // sequential
  for(i=0;i<4;i++)hmvm_cuda0<double>(matD2, b, i, dump_result, nbench);
  // block
  for(i=0;i<4;i++)hmvm_cuda0<double>(matD2, b, 10+i, dump_result, nbench);
  // hybrid1: DIV(1,2,4,8,16,32), ATOMIC(1,2), = 96patterns
  for(i=0;i<96;i++)hmvm_cuda1<double>(matD2, b, i, dump_result, nbench);
  // hybrid2: DIV(1,2,4,8,16,32), MUL(1,2,3,...,16), ATOMIC(1,2)=1536patterns
  for(i=0;i<1536;i++)hmvm_cuda2<double>(matD2, b, i, dump_result, nbench);
  // hybrid3: DIV(1,2,4,8,16,32), MUL(1,2,3,...,16), ATOMIC(1,2)=1536patterns
  for(i=0;i<1536;i++)hmvm_cuda3<double>(matD2, b, i, dump_result, nbench);
#endif
  // MAGMA BLAS
  hmvm_magma<double>(matD2, b, 0, dump_result, nbench);
  hmvm_magma_batched1<double>(matD2, b, 0, dump_result, nbench);
  hmvm_magma_batched1<double>(matD2, b, 1, dump_result, nbench);
  hmvm_magma_batched2<double>(matD2, b, 0, dump_result, nbench);
  hmvm_magma_batched2<double>(matD2, b, 1, dump_result, nbench);
  delete [] b;
#endif

#if 1
  // ######## ######## ######## ######## ######## ######## ######## ########
  // float
  // ######## ######## ######## ######## ######## ######## ######## ########
  printf("######## ######## ######## ########\n");
  printf("float\n");
  printf("######## ######## ######## ########\n");
  matF = new matrix<float>();
  matF2 = new matrix2<float>();
  float *fb=NULL;
  fb = new float[nd];
  for(i=0;i<nd;i++){
	fb[i] = sin((float)(i+1));
  }
  convertD2F(matF, matF2, matD, matD2);
  // sequential
  hmvm_seq<float>(matF, matF2, fb, dump_result, nbench);
  // OpenMP
  hmvm_omp<float>(matF, matF2, fb, dump_result, nbench);
  // CUDA
#ifdef _USE_CUDA
  // sequential
  for(i=0;i<4;i++)hmvm_cuda0<float>(matF2, fb, i, dump_result, nbench);
  // block
  for(i=0;i<4;i++)hmvm_cuda0<float>(matF2, fb, 10+i, dump_result, nbench);
  // hybrid1: DIV(1,2,4,8,16,32), ATOMIC(1,2), = 96patterns
  for(i=0;i<96;i++)hmvm_cuda1<float>(matF2, fb, i, dump_result, nbench);
  // hybrid2: DIV(1,2,4,8,16,32), MUL(1,2,3,...,16), ATOMIC(1,2)=1536patterns
  for(i=0;i<1536;i++)hmvm_cuda2<float>(matF2, fb, i, dump_result, nbench);
  // hybrid3: DIV(1,2,4,8,16,32), MUL(1,2,3,...,16), ATOMIC(1,2)=1536patterns
  for(i=0;i<1536;i++)hmvm_cuda3<float>(matF2, fb, i, dump_result, nbench);
#endif
  // MAGMA BLAS
  // MAGMA BLAS
  hmvm_magma<float>(matF2, fb, 0, dump_result, nbench);
  hmvm_magma_batched<float>(matF2, fb, 0, dump_result, nbench);
  hmvm_magma_batched<float>(matF2, fb, 1, dump_result, nbench);
  hmvm_magma_batched2<float>(matF2, fb, 0, dump_result, nbench);
  hmvm_magma_batched2<float>(matF2, fb, 1, dump_result, nbench);
  delete [] fb;
#endif

  delete [] matF2;
  delete [] matF;
  delete [] matD2;
  delete [] matD;
  return 0;
}
