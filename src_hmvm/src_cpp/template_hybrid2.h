#ifndef _TEMPLATE_HYBRID2_H
#define _TEMPLATE_HYBRID2_H

#if 0

template __global__ void hmvm_cuda_hybrid2<float,DIV,MUL,ATOMIC>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
int napprox, int *approx, int ndense, int *dense);

for a in 0 1
do
for m in `seq 1 16`
do
for d in 1 2 4 8 16 32
do
echo "template __global__ void hmvm_cuda_hybrid2<float,${d},${m},${a}>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
int napprox, int *approx, int ndense, int *dense);"
done
done
done
*/

/*
 template __global__ void hmvm_cuda_hybrid2<double,DIV,MUL,ATOMIC>
(T *d_zaut, T *d_zu, int nlf, int ktmax,
int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
int napprox, int *approx, int ndense, int *dense);

for a in 0 1
do
for m in `seq 1 16`
do
for d in 1 2 4 8 16 32
do
echo "template __global__ void hmvm_cuda_hybrid2<double,${d},${m},${a}>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
int napprox, int *approx, int ndense, int *dense);"
done
done
done


#endif
