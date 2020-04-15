#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hacapk.h"

int main(int argc, char **argv)
{
  int i, ret;
  FILE *F;
  matrix *mat=NULL;
  matrix2 *mat2=NULL;
  int dump_result = 0;

  if(argc<2){
	printf("usage: %s filename [dump_result]\n", argv[0]);
	return -1;
  }
  if(argc>=3){dump_result=atoi(argv[2]);printf("dump_result = %d\n",dump_result);}
  printf("call loadHmatrix %s\n", argv[1]); fflush(stdout);
  mat = (matrix*)malloc(sizeof(matrix));
  mat2 = (matrix2*)malloc(sizeof(matrix2));
  ret = loadHmatrix(argv[1], mat, mat2);
  if(ret)return -1;

  int nd = mat->nd;
  printf("nd=%d\n",nd);
  double *b=NULL, *v=NULL;
  b=(double*)malloc(sizeof(double)*nd);
  v=(double*)malloc(sizeof(double)*nd);
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
  return 0;
}
