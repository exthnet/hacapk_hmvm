#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hacapk.h"
#include "hmvm_seq.h"
#include "hmvm_omp.h"
#include "hmvm_cuda.h"

int main(int argc, char **argv)
{
  int i, ret;
  FILE *F;
  matrix<double> *matD=NULL;
  matrix2<double> *matD2=NULL;
  matrix<float> *matF=NULL;
  matrix2<float> *matF2=NULL;
  int dump_result=0;

  if(argc<2){
	printf("usage: %s filename [opt]\n", argv[0]);
	return -1;
  }
  if(argc>=3){dump_result=atoi(argv[2]);printf("dump_result = %d\n",dump_result);}
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
  hmvm_seq<double>(matD, matD2, b, dump_result);
  hmvm_seq_bench<double>(matD, matD2, b);
  // OpenMP
  hmvm_omp<double>(matD, matD2, b, dump_result);
  hmvm_omp_bench<double>(matD, matD2, b);
  // MKL
#ifdef __INTEL_COMPILER
  hmvm_blas_p(matD, matD2, b, dump_result);
  hmvm_blas_p_bench(matD, matD2, b);
  hmvm_blas_s(matD, matD2, b, dump_result);
  hmvm_blas_s_bench(matD, matD2, b);
  hmvm_cblas_p(matD, matD2, b, dump_result);
  hmvm_cblas_p_bench(matD, matD2, b);
  hmvm_cblas_s(matD, matD2, b, dump_result);
  hmvm_cblas_s_bench(matD, matD2, b);
#endif
  // CUDA
#ifdef _USE_CUDA
  hmvm_cuda1<double>(matD2, b, 0, dump_result);
  hmvm_cuda1<double>(matD2, b, 1, dump_result);
  hmvm_cuda1<double>(matD2, b, 2, dump_result);
  hmvm_cuda1<double>(matD2, b, 10, dump_result); // 1
  hmvm_cuda1<double>(matD2, b, 11, dump_result); // 2
  hmvm_cuda1<double>(matD2, b, 12, dump_result); // 4
  hmvm_cuda1<double>(matD2, b, 13, dump_result); // 8
  hmvm_cuda1<double>(matD2, b, 14, dump_result); // 16
  hmvm_cuda1<double>(matD2, b, 15, dump_result); // 32
  hmvm_cuda1<double>(matD2, b, 20, dump_result); // 1
  hmvm_cuda1<double>(matD2, b, 21, dump_result); // 2
  hmvm_cuda1<double>(matD2, b, 22, dump_result); // 4
  hmvm_cuda1<double>(matD2, b, 23, dump_result); // 8
  hmvm_cuda1<double>(matD2, b, 24, dump_result); // 16
  hmvm_cuda1<double>(matD2, b, 25, dump_result); // 32
#endif
  delete [] b;//free(b);
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
  fb = new float[nd];//(float*)malloc(sizeof(float)*matD.nd);
  for(i=0;i<nd;i++){
	fb[i] = sin((float)(i+1));
  }
  convertD2F(matF, matF2, matD, matD2);
  // sequential
  hmvm_seq<float>(matF, matF2, fb, dump_result);
  hmvm_seq_bench<float>(matF, matF2, fb);
  // OpenMP
  hmvm_omp<float>(matF, matF2, fb, dump_result);
  hmvm_omp_bench<float>(matF, matF2, fb);
  // CUDA
#ifdef _USE_CUDA
  hmvm_cuda1<float>(matF2, fb, 0, dump_result);
  hmvm_cuda1<float>(matF2, fb, 1, dump_result);
  hmvm_cuda1<float>(matF2, fb, 2, dump_result);
  hmvm_cuda1<float>(matF2, fb, 10, dump_result); // 1
  hmvm_cuda1<float>(matF2, fb, 11, dump_result); // 2
  hmvm_cuda1<float>(matF2, fb, 12, dump_result); // 4
  hmvm_cuda1<float>(matF2, fb, 13, dump_result); // 8
  hmvm_cuda1<float>(matF2, fb, 14, dump_result); // 16
  hmvm_cuda1<float>(matF2, fb, 15, dump_result); // 32
  hmvm_cuda1<float>(matF2, fb, 20, dump_result); // 1
  hmvm_cuda1<float>(matF2, fb, 21, dump_result); // 2
  hmvm_cuda1<float>(matF2, fb, 22, dump_result); // 4
  hmvm_cuda1<float>(matF2, fb, 23, dump_result); // 8
  hmvm_cuda1<float>(matF2, fb, 24, dump_result); // 16
  hmvm_cuda1<float>(matF2, fb, 25, dump_result); // 32
#endif
  delete [] fb;//free(fb);
#endif
  delete [] matF2;
  delete [] matF;
  delete [] matD2;
  delete [] matD;
  return 0;
}
