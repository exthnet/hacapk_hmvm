#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hacapk.h"

int loadHmatrix(const char *fname, matrix<double> *mat, matrix2<double> *mat2)
{
  int i, ret;
  FILE *F;
  int irecord;
  double drecord;
  int nd, nlf, ktmax, len, offset;
  printf("loadHmatrix: begin\n"); fflush(stdout);
  F = fopen(fname, "r");
  if(F==NULL){
	printf("fopen %s failed\n", fname);
	return -1;
  }
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(&nd, sizeof(int), 1, F);
  ret = fread(&nlf, sizeof(int), 1, F);
  ret = fread(&ktmax, sizeof(int), 1, F);
  ret = fread(&irecord, sizeof(int), 1, F);
  printf("nd =  %d\n", nd);
  printf("nlf = %d\n", nlf);
  printf("ktmax = %d\n", ktmax);
  mat->nd = nd;
  mat->nlf = nlf;
  mat->ktmax = ktmax;

  mat->submat = (submatrix<double>*)malloc(sizeof(submatrix<double>)*mat->nlf);
  if(mat->submat==NULL){
	printf("malloc mat->submat failed\n");
	return -1;
  }

  len = 0;

  printf("loadHmatrix: load rowdata\n"); fflush(stdout);
  for(i=0;i<mat->nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	ret = fread(&irecord, sizeof(int), 1, F);
	ret = fread(&ltmtx, sizeof(int), 1, F);
	ret = fread(&ndl, sizeof(int), 1, F);
	ret = fread(&ndt, sizeof(int), 1, F);
	ret = fread(&nstrtl, sizeof(int), 1, F);
	ret = fread(&nstrtt, sizeof(int), 1, F);
	ret = fread(&kt, sizeof(int), 1, F);
	ret = fread(&irecord, sizeof(int), 1, F);
	//printf("%d: %d %d %d %d %d %d\n", i,ltmtx, ndl, ndt, nstrtl, nstrtt, kt);

	mat->submat[i].ltmtx = ltmtx;
	mat->submat[i].ndl = ndl;
	mat->submat[i].ndt = ndt;
	mat->submat[i].nstrtl = nstrtl;
	mat->submat[i].nstrtt = nstrtt;
	mat->submat[i].kt = kt;

	if(ltmtx==1){
	  mat->submat[i].a1 = (double*)malloc(sizeof(double)*kt*ndt);
	  mat->submat[i].a2 = (double*)malloc(sizeof(double)*kt*ndl);
	  mat->submat[i].a2t = (double*)malloc(sizeof(double)*kt*ndl);
	  ret = fread(&irecord, sizeof(int), 1, F);
	  ret = fread(mat->submat[i].a1, sizeof(double), kt*ndt, F);
	  ret = fread(&irecord, sizeof(int), 1, F);
	  ret = fread(&irecord, sizeof(int), 1, F);
	  ret = fread(mat->submat[i].a2, sizeof(double), kt*ndl, F);
	  ret = fread(&irecord, sizeof(int), 1, F);
	  len += kt*ndt + kt*ndl;
	  int x, y;
	  for(y=0;y<kt;y++){
		for(x=0;x<ndl;x++){
		  mat->submat[i].a2t[y+x*kt] = mat->submat[i].a2[y*ndl+x];
		}
	  }
	}else{
	  mat->submat[i].a1 = (double*)malloc(sizeof(double)*ndl*ndt);
	  mat->submat[i].a2 = NULL;
	  mat->submat[i].a2t = NULL;
	  ret = fread(&irecord, sizeof(int), 1, F);
	  ret = fread(mat->submat[i].a1, sizeof(double), ndl*ndt, F);
	  ret = fread(&irecord, sizeof(int), 1, F);
	  len += ndl*ndt;
	}
  }

  printf("loadHmatrix: make matrix2\n"); fflush(stdout);
  mat2->nd = nd;
  mat2->nlf = nlf;
  mat2->ktmax = ktmax;
  mat2->ltmtx = (int*)malloc(sizeof(int)*nlf);
  mat2->kt = (int*)malloc(sizeof(int)*nlf);
  mat2->ndl = (int*)malloc(sizeof(int)*nlf);
  mat2->ndt = (int*)malloc(sizeof(int)*nlf);
  mat2->nstrtl = (int*)malloc(sizeof(int)*nlf);
  mat2->nstrtt = (int*)malloc(sizeof(int)*nlf);
  mat2->a1 = (int*)malloc(sizeof(int)*nlf);
  mat2->a2 = (int*)malloc(sizeof(int)*nlf);
  mat2->rowmat = (double*)malloc(sizeof(double)*len);
  mat2->rowmat_t = (double*)malloc(sizeof(double)*len);
  offset = 0;
  for(i=0;i<nlf;i++){
	mat2->ltmtx[i] = mat->submat[i].ltmtx;
	mat2->kt[i] = mat->submat[i].kt;
	mat2->ndl[i] = mat->submat[i].ndl;
	mat2->ndt[i] = mat->submat[i].ndt;
	mat2->nstrtl[i] = mat->submat[i].nstrtl;
	mat2->nstrtt[i] = mat->submat[i].nstrtt;
	int ltmtx = mat->submat[i].ltmtx;
	if(ltmtx==1){
	  mat2->a1[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndt);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndt);
	  offset += mat->submat[i].kt*mat->submat[i].ndt;
	  mat2->a2[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a2, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndl);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a2t, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndl);
	  offset += mat->submat[i].kt*mat->submat[i].ndl;
	}else{
	  mat2->a1[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].ndl*mat->submat[i].ndt);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].ndl*mat->submat[i].ndt);
	  offset += mat->submat[i].ndl*mat->submat[i].ndt;
	  mat2->a2[i] = 0;
	}
  }
  mat2->len = offset;

  printf("loadHmatrix: end\n"); fflush(stdout);

  // for debug
#if 0
  for(i=0;i<mat->nlf;i++){
	int j;
	if(mat->submat[i].ltmtx==0){
	  for(j=0; j<mat->submat[i].kt*mat->submat[i].ndt; j++)mat->submat[i].a1[j] = 1.0;
	  for(j=0; j<mat->submat[i].kt*mat->submat[i].ndl; j++)mat->submat[i].a2[j] = 1.0;
	  for(j=0; j<mat->submat[i].kt*mat->submat[i].ndl; j++)mat->submat[i].a2t[j] = 1.0;
	}else{
	  for(j=0; j<mat->submat[i].ndl*mat->submat[i].ndt; j++)mat->submat[i].a1[j] = 1.0;
	}
  }
  for(i=0; i<len; i++){
	mat2->rowmat[i] = 1.0;
	mat2->rowmat_t[i] = 1.0;
  }
#endif
  return 0;
}

int loadHmatrix2(const char *fname, matrix2<double> *mat2)
{
  int i, ret;
  FILE *F;
  int irecord;
  int nd, nlf, ktmax, len, offset;
  printf("loadHmatrix: begin\n"); fflush(stdout);
  F = fopen(fname, "r");
  if(F==NULL){
	printf("fopen %s failed\n", fname);
	return -1;
  }
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(&nd, sizeof(int), 1, F);
  ret = fread(&nlf, sizeof(int), 1, F);
  ret = fread(&ktmax, sizeof(int), 1, F);
  ret = fread(&len, sizeof(int), 1, F);
  ret = fread(&irecord, sizeof(int), 1, F);
  printf("nd    = %d\n", nd);
  printf("nlf   = %d\n", nlf);
  printf("ktmax = %d\n", ktmax);
  printf("len   = %d\n", len);
  fflush(stdout);
  mat2->nd = nd;
  mat2->nlf = nlf;
  mat2->ktmax = ktmax;
  mat2->len = len;

  printf("loadHmatrix: make matrix2\n"); fflush(stdout);
  mat2->ltmtx = (int*)malloc(sizeof(int)*nlf);
  mat2->kt = (int*)malloc(sizeof(int)*nlf);
  mat2->ndl = (int*)malloc(sizeof(int)*nlf);
  mat2->ndt = (int*)malloc(sizeof(int)*nlf);
  mat2->nstrtl = (int*)malloc(sizeof(int)*nlf);
  mat2->nstrtt = (int*)malloc(sizeof(int)*nlf);
  mat2->a1 = (int*)malloc(sizeof(int)*nlf);
  mat2->a2 = (int*)malloc(sizeof(int)*nlf);

  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(mat2->ltmtx, sizeof(int), nlf, F); printf("ltmtx %d\n", mat2->ltmtx[0]); fflush(stdout);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(mat2->ndl, sizeof(int), nlf, F); printf("ndl %d\n", mat2->ndl[0]); fflush(stdout);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(mat2->ndt, sizeof(int), nlf, F); printf("ndt %d\n", mat2->ndt[0]); fflush(stdout);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(mat2->nstrtl, sizeof(int), nlf, F); printf("nstrtl %d\n", mat2->nstrtl[0]); fflush(stdout);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(mat2->nstrtt, sizeof(int), nlf, F); printf("nstrtt %d\n", mat2->nstrtt[0]); fflush(stdout);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(&irecord, sizeof(int), 1, F);
  ret = fread(mat2->kt, sizeof(int), nlf, F); printf("kt %d\n", mat2->kt[0]); fflush(stdout);
  ret = fread(&irecord, sizeof(int), 1, F);

  mat2->rowmat = (double*)malloc(sizeof(double)*len);
  mat2->rowmat_t = (double*)malloc(sizeof(double)*len);

  offset = 0;
  for(i=0;i<nlf;i++){
	int ltmtx = mat2->ltmtx[i];
	if(ltmtx==1){
	  mat2->a1[i] = offset;
	  ret = fread(&irecord, sizeof(int), 1, F);
	  ret = fread(&mat2->rowmat[offset], sizeof(double), mat2->kt[i]*mat2->ndt[i], F);
	  ret = fread(&irecord, sizeof(int), 1, F);
	  memcpy(&mat2->rowmat_t[offset], &mat2->rowmat[offset], sizeof(double)*mat2->kt[i]*mat2->ndt[i]);
	  offset += mat2->kt[i]*mat2->ndt[i];
	  mat2->a2[i] = offset;
	  ret = fread(&irecord, sizeof(int), 1, F);
	  ret = fread(&mat2->rowmat[offset], sizeof(double), mat2->kt[i]*mat2->ndl[i], F);
	  ret = fread(&irecord, sizeof(int), 1, F);
	  //memcpy(&mat2->rowmat_t[offset], mat->submat[i].a2t, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndl);
	  int x, y;
	  int kt = mat2->kt[i];
	  int ndl = mat2->ndl[i];
	  for(y=0;y<kt;y++){
		for(x=0;x<ndl;x++){
		  mat2->rowmat_t[offset+y+x*kt] = mat2->rowmat[offset+y*ndl+x];
		}
	  }
	  offset += mat2->kt[i]*mat2->ndl[i];
	}else{
	  mat2->a1[i] = offset;
	  ret = fread(&irecord, sizeof(int), 1, F);
	  ret = fread(&mat2->rowmat[offset], sizeof(double), mat2->ndl[i]*mat2->ndt[i], F);
	  ret = fread(&irecord, sizeof(int), 1, F);
	  memcpy(&mat2->rowmat_t[offset], &mat2->rowmat[offset], sizeof(double)*mat2->ndl[i]*mat2->ndt[i]);
	  offset += mat2->ndl[i]*mat2->ndt[i];
	  mat2->a2[i] = 0;
	}
  }
  printf("check len: %d %d\n", mat2->len, offset);

  printf("loadHmatrix2: end\n"); fflush(stdout);
  return 0;
}

int convertD2F(matrix<float> *matF, matrix2<float> *matF2, const matrix<double> *matD, const matrix2<double> *matD2)
{
  int i, j;
  printf("convertD2F: begin\n"); fflush(stdout);
  if(matD!=NULL){
	int len;
	printf("convertD2F: matD -> matF: begin\n"); fflush(stdout);
	matF->nd     = matD->nd;
	matF->nlf    = matD->nlf;
	matF->ktmax  = matD->ktmax;
	matF->submat = new submatrix<float>[matF->nlf]; //(submatrix<float>*)malloc(sizeof(submatrix<float>)*matF->nlf);
	len = 0;
	for(i=0;i<matF->nlf;i++){
	  matF->submat[i].ltmtx = matD->submat[i].ltmtx;
	  matF->submat[i].ndl    = matD->submat[i].ndl;
	  matF->submat[i].ndt    = matD->submat[i].ndt;
	  matF->submat[i].nstrtl = matD->submat[i].nstrtl;
	  matF->submat[i].nstrtt = matD->submat[i].nstrtt;
	  matF->submat[i].kt     = matD->submat[i].kt;
	  int kt  = matD->submat[i].kt;
	  int ndt = matD->submat[i].ndt;
	  int ndl = matD->submat[i].ndl;
	  if(matF->submat[i].ltmtx==1){
		matF->submat[i].a1  = new float[kt*ndt]; //(float*)malloc(sizeof(float)*kt*ndt);
		matF->submat[i].a2  = new float[kt*ndl]; //(float*)malloc(sizeof(float)*kt*ndl);
		matF->submat[i].a2t = new float[kt*ndl]; //(float*)malloc(sizeof(float)*kt*ndl);
		for(j=0; j<kt*ndt; j++)matF->submat[i].a1[j]  = matD->submat[i].a1[j];
		for(j=0; j<kt*ndl; j++)matF->submat[i].a2[j]  = matD->submat[i].a2[j];
		for(j=0; j<kt*ndl; j++)matF->submat[i].a2t[j] = matD->submat[i].a2t[j];
		len += kt*ndt + kt*ndl;
	  }else{
		matF->submat[i].a1  = new float[ndl*ndt]; //(float*)malloc(sizeof(float)*ndl*ndt);
		matF->submat[i].a2  = NULL;
		matF->submat[i].a2t = NULL;
		for(j=0; j<ndl*ndt; j++)matF->submat[i].a1[j] = matD->submat[i].a1[j];
		len += ndl*ndt;
	  }
	}
	printf("convertD2F: matD -> matF: end\n"); fflush(stdout);
  }
  if(matD2!=NULL){
	int offset;
	printf("convertD2F: matD2 -> matF2: begin\n"); fflush(stdout);
	matF2->nd       = matD2->nd;
	matF2->nlf      = matD2->nlf;
	matF2->ktmax    = matD2->ktmax;
	matF2->ltmtx    = new int[matD2->nlf]; //(int*)malloc(sizeof(int)*matD2->nlf);
	matF2->kt       = new int[matD2->nlf]; //(int*)malloc(sizeof(int)*matD2->nlf);
	matF2->ndl      = new int[matD2->nlf]; //(int*)malloc(sizeof(int)*matD2->nlf);
	matF2->ndt      = new int[matD2->nlf]; //(int*)malloc(sizeof(int)*matD2->nlf);
	matF2->nstrtl   = new int[matD2->nlf]; //(int*)malloc(sizeof(int)*matD2->nlf);
	matF2->nstrtt   = new int[matD2->nlf]; //(int*)malloc(sizeof(int)*matD2->nlf);
	matF2->a1       = new int[matD2->nlf]; //(int*)malloc(sizeof(int)*matD2->nlf);
	matF2->a2       = new int[matD2->nlf]; //(int*)malloc(sizeof(int)*matD2->nlf);
	matF2->rowmat   = new float[matD2->len]; //(float*)malloc(sizeof(float)*len);
	matF2->rowmat_t = new float[matD2->len]; //(float*)malloc(sizeof(float)*len);
	offset = 0;
	for(i=0;i<matD2->nlf;i++){
	  matF2->ltmtx[i]  = matD2->ltmtx[i];
	  matF2->ndl[i]    = matD2->ndl[i];
	  matF2->ndt[i]    = matD2->ndt[i];
	  matF2->nstrtl[i] = matD2->nstrtl[i];
	  matF2->nstrtt[i] = matD2->nstrtt[i];
	  matF2->kt[i]     = matD2->kt[i];
	  int kt  = matD2->kt[i];
	  int ndt = matD2->ndt[i];
	  int ndl = matD2->ndl[i];
	  if(matF2->ltmtx[i]==1){
		matF2->a1[i] = offset;
		for(j=0;j<kt*ndt;j++)matF2->rowmat[offset+j]   = matD2->rowmat[offset+j];
		for(j=0;j<kt*ndt;j++)matF2->rowmat_t[offset+j] = matD2->rowmat_t[offset+j];
		offset += kt*ndt;
		matF2->a2[i] = offset;
		for(j=0;j<kt*ndl;j++)matF2->rowmat[offset+j]   = matD2->rowmat[offset+j];
		for(j=0;j<kt*ndl;j++)matF2->rowmat_t[offset+j] = matD2->rowmat_t[offset+j];
		offset += kt*ndl;
	  }else{
		matF2->a1[i] = offset;
		for(j=0;j<ndl*ndt;j++)matF2->rowmat[offset+j]   = matD2->rowmat[offset+j];
		for(j=0;j<ndl*ndt;j++)matF2->rowmat_t[offset+j] = matD2->rowmat_t[offset+j];
		offset += ndl*ndt;
		matF2->a2[i] = 0;
	  }
	}
	matF2->len = matD2->len;
	printf("check len: %d %d\n", matF2->len, offset);
	printf("convertD2F: matD2 -> matF2: end\n"); fflush(stdout);
  }
  printf("convertD2F: end\n"); fflush(stdout);
  return 0;
}

