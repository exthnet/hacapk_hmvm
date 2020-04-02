#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hacapk_c.h"

int main(int argc, char **argv)
{
  int i, ret;
  FILE *F;
  matrix mat;
  matrix2 mat2;
  int dump_result = 0;

  if(argc<2){
	printf("usage: %s filename [dump_result]\n", argv[0]);
	return -1;
  }
  if(argc>=3){dump_result=atoi(argv[2]);printf("dump_result = %d\n",dump_result);}
  printf("call loadmatrix %s\n", argv[1]); fflush(stdout);
  ret = loadmatrix(argv[1], &mat, &mat2);
  if(ret)return -1;

  int nd = mat.nd;
  printf("nd=%d\n",nd);
  double *b=NULL, *v=NULL;
  b=(double*)malloc(sizeof(double)*mat.nd);
  v=(double*)malloc(sizeof(double)*mat.nd);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
	v[i] = 0.0;
  }

  F = fopen("vec.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", b[i]);
  fclose(F);

  // ######## ######## ######## ######## ######## ######## ######## ########
  // ready for hmvm
  // ######## ######## ######## ######## ######## ######## ######## ########
#if 0
  int ip;
  mat2.approx = (int*)malloc(sizeof(int)*mat2.nlf);
  mat2.dense  = (int*)malloc(sizeof(int)*mat2.nlf);
  mat2.napprox = mat2.ndense = 0;
  for(ip=0; ip<mat2.nlf; ip++){
	if(mat2.ltmtx[ip]==1){
	  mat2.approx[mat2.napprox++] = ip;
	}else{
	  mat2.dense[mat2.ndense++] = ip;
	}
  }
  printf(" %d approx, %d dense\n", mat2.napprox, mat2.ndense);
#endif
  // sequential
  hmvm_seq(mat, mat2, b, dump_result);
  hmvm_seq_bench(mat, mat2, b);
  // OpenMP
  hmvm_omp(mat, mat2, b, dump_result);
  hmvm_omp_bench(mat, mat2, b);
  // MKL
#ifdef __INTEL_COMPILER
  hmvm_blas_p(mat, mat2, b, dump_result);
  hmvm_blas_p_bench(mat, mat2, b);
  hmvm_blas_s(mat, mat2, b, dump_result);
  hmvm_blas_s_bench(mat, mat2, b);
  hmvm_cblas_p(mat, mat2, b, dump_result);
  hmvm_cblas_p_bench(mat, mat2, b);
  hmvm_cblas_s(mat, mat2, b, dump_result);
  hmvm_cblas_s_bench(mat, mat2, b);
#endif
  // CUDA
#ifdef _USE_CUDA
  hmvm_cuda1(mat2, b, 0, dump_result);
  hmvm_cuda1(mat2, b, 1, dump_result);
#endif
  return 0;
}
