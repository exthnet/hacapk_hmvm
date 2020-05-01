#ifndef _TEMPLATE_HYBRID1_H
#define _TEMPLATE_HYBRID1_H

template void hmvm_cuda_hybrid1_proxy<float,auto,auto,auto,auto,auto>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, double *v, double *b, int nd, char *fname, int bench);

#endif
