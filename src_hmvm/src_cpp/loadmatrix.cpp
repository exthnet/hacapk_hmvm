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
  return 0;
}

int convertD2F(matrix<float> *matF, matrix2<float> *matF2, matrix<double> matD, matrix2<double> matD2)
{
  int i, j;
  int len, offset;
  matF->nd     = matD.nd;
  matF->nlf    = matD.nlf;
  matF->ktmax  = matD.ktmax;
  matF->submat = new submatrix<float>[matF->nlf]; //(submatrix<float>*)malloc(sizeof(submatrix<float>)*matF->nlf);
  len = 0;
  for(i=0;i<matF->nlf;i++){
	matF->submat[i].ltmtx = matD.submat[i].ltmtx;
	matF->submat[i].ndl    = matD.submat[i].ndl;
	matF->submat[i].ndt    = matD.submat[i].ndt;
	matF->submat[i].nstrtl = matD.submat[i].nstrtl;
	matF->submat[i].nstrtt = matD.submat[i].nstrtt;
	matF->submat[i].kt     = matD.submat[i].kt;
	int kt  = matD.submat[i].kt;
	int ndt = matD.submat[i].ndt;
	int ndl = matD.submat[i].ndl;
	if(matF->submat[i].ltmtx==1){
	  matF->submat[i].a1  = new float[kt*ndt]; //(float*)malloc(sizeof(float)*kt*ndt);
	  matF->submat[i].a2  = new float[kt*ndl]; //(float*)malloc(sizeof(float)*kt*ndl);
	  matF->submat[i].a2t = new float[kt*ndl]; //(float*)malloc(sizeof(float)*kt*ndl);
	  for(j=0; j<kt*ndt; j++)matF->submat[i].a1[j]  = matD.submat[i].a1[j];
	  for(j=0; j<kt*ndl; j++)matF->submat[i].a2[j]  = matD.submat[i].a2[j];
	  for(j=0; j<kt*ndl; j++)matF->submat[i].a2t[j] = matD.submat[i].a2t[j];
	  len += kt*ndt + kt*ndl;
	}else{
	  matF->submat[i].a1  = new float[ndl*ndt]; //(float*)malloc(sizeof(float)*ndl*ndt);
	  matF->submat[i].a2  = NULL;
	  matF->submat[i].a2t = NULL;
	  for(j=0; j<ndl*ndt; j++)matF->submat[i].a1[j] = matD.submat[i].a1[j];
	  len += ndl*ndt;
	}
  }

  matF2->nd       = matD2.nd;
  matF2->nlf      = matD2.nlf;
  matF2->ktmax    = matD2.ktmax;
  matF2->ltmtx    = new int[matD2.nlf]; //(int*)malloc(sizeof(int)*matD2.nlf);
  matF2->kt       = new int[matD2.nlf]; //(int*)malloc(sizeof(int)*matD2.nlf);
  matF2->ndl      = new int[matD2.nlf]; //(int*)malloc(sizeof(int)*matD2.nlf);
  matF2->ndt      = new int[matD2.nlf]; //(int*)malloc(sizeof(int)*matD2.nlf);
  matF2->nstrtl   = new int[matD2.nlf]; //(int*)malloc(sizeof(int)*matD2.nlf);
  matF2->nstrtt   = new int[matD2.nlf]; //(int*)malloc(sizeof(int)*matD2.nlf);
  matF2->a1       = new int[matD2.nlf]; //(int*)malloc(sizeof(int)*matD2.nlf);
  matF2->a2       = new int[matD2.nlf]; //(int*)malloc(sizeof(int)*matD2.nlf);
  matF2->rowmat   = new float[len]; //(float*)malloc(sizeof(float)*len);
  matF2->rowmat_t = new float[len]; //(float*)malloc(sizeof(float)*len);
  offset = 0;
  for(i=0;i<matD2.nlf;i++){
	matF2->ltmtx[i]  = matF->submat[i].ltmtx;
	matF2->ndl[i]    = matF->submat[i].ndl;
	matF2->ndt[i]    = matF->submat[i].ndt;
	matF2->nstrtl[i] = matF->submat[i].nstrtl;
	matF2->nstrtt[i] = matF->submat[i].nstrtt;
	matF2->kt[i]     = matF->submat[i].kt;
	int kt  = matF->submat[i].kt;
	int ndt = matF->submat[i].ndt;
	int ndl = matF->submat[i].ndl;
	if(matF2->ltmtx[i]==1){
	  matF2->a1[i] = offset;
	  memcpy(&matF2->rowmat[offset], matF->submat[i].a1, sizeof(float)*kt*ndt);
	  memcpy(&matF2->rowmat_t[offset], matF->submat[i].a1, sizeof(float)*kt*ndt);
	  offset += kt*ndt;
	  matF2->a2[i] = offset;
	  memcpy(&matF2->rowmat[offset], matF->submat[i].a2, sizeof(float)*kt*ndl);
	  memcpy(&matF2->rowmat_t[offset], matF->submat[i].a2t, sizeof(float)*kt*ndl);
	  offset += kt*ndl;
	}else{
	  matF2->a1[i] = offset;
	  memcpy(&matF2->rowmat[offset], matF->submat[i].a1, sizeof(double)*ndl*ndt);
	  memcpy(&matF2->rowmat_t[offset], matF->submat[i].a1, sizeof(double)*ndl*ndt);
	  offset += ndl*ndt;
	  matF2->a2[i] = 0;
	}
  }
  matF2->len = offset;

  return 0;
}

