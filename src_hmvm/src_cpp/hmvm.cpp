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
  matrix<double> matD;
  matrix2<double> matD2;
  matrix<float> matF;
  matrix2<float> matF2;
  int dump_result=0;

  if(argc<2){
	printf("usage: %s filename [opt]\n", argv[0]);
	return -1;
  }
  if(argc>=3){dump_result=atoi(argv[2]);printf("dump_result = %d\n",dump_result);}
  printf("call loadmatrix %s\n", argv[1]); fflush(stdout);
  ret = loadHmatrix(argv[1], &matD, &matD2);
  if(ret)return -1;

  int nd = matD.nd;
  printf("nd=%d\n",nd);
  double *b=NULL;
  b = new double[matD.nd];//(double*)malloc(sizeof(double)*matD.nd);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
  }

  F = fopen("vec.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", b[i]);
  fclose(F);

  // ######## ######## ######## ######## ######## ######## ######## ########
  // ready for hmvm
  // ######## ######## ######## ######## ######## ######## ######## ########
#if 0
  int ip;
  matD2.approx = (int*)malloc(sizeof(int)*matD2.nlf);
  matD2.dense  = (int*)malloc(sizeof(int)*matD2.nlf);
  matD2.napprox = matD2.ndense = 0;
  for(ip=0; ip<matD2.nlf; ip++){
	if(matD2.ltmtx[ip]==1){
	  matD2.approx[matD2.napprox++] = ip;
	}else{
	  matD2.dense[matD2.ndense++] = ip;
	}
  }
  printf(" %d approx, %d dense\n", matD2.napprox, matD2.ndense);
#endif

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
  hmvm_cuda1(matD2, b, 0, dump_result);
  hmvm_cuda1(matD2, b, 1, dump_result);
  hmvm_cuda1(matD2, b, 2, dump_result);
#endif
  delete [] b;//free(b);

  // ######## ######## ######## ######## ######## ######## ######## ########
  // float
  // ######## ######## ######## ######## ######## ######## ######## ########
  printf("######## ######## ######## ########\n");
  printf("float\n");
  printf("######## ######## ######## ########\n");
  float *fb=NULL;
  fb = new float[matD.nd];//(float*)malloc(sizeof(float)*matD.nd);
  for(i=0;i<nd;i++){
	fb[i] = sin((float)(i+1));
  }
  convertD2F(&matF, &matF2, matD, matD2);
  // sequential
  hmvm_seq<float>(matF, matF2, fb, dump_result);
  hmvm_seq_bench<float>(matF, matF2, fb);
  // OpenMP
  hmvm_omp<float>(matF, matF2, fb, dump_result);
  hmvm_omp_bench<float>(matF, matF2, fb);
  // CUDA
#ifdef _USE_CUDA
  hmvm_cuda1(matF2, fb, 0, dump_result);
  hmvm_cuda1(matF2, fb, 1, dump_result);
  hmvm_cuda1(matF2, fb, 2, dump_result);
#endif
  delete [] fb;//free(fb);

  return 0;
}
