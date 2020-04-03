# hacapk_hmvm

## src_hdump
### 概要
HACApK 1.0.0を改造し、H行列をバイナリファイルに出力するようにしたもの。
m_HACApK_hdump.f90以外のコードの多くはHACApK 1.0.0のコードを流用している。
### 使い方
makeでbem-bb-SCM.outファイルが生成される。
MPI + C 環境であればどこでも動く、と期待される。

./bem-bb-SCM.out inputmatrix
として実行すると、hmatrix_0.binというバイナリファイルが出力される。
出力されたバイナリファイルはsrc_hmvmの入力ファイルとして使える。

## src_hmvm
### 概要
src_hdumpで出力したバイナリファイルを用いてH行列ベクトル積を行うテストコード群。
### 使い方

- src_c
シンプルなC言語コード、CPU(OpenMP含む)向け、倍精度。

- src_cpp
C++コード、CPU(OpenMP)とGPU(CUDA)、単精度と倍精度。
実装(移植)中。
-- Makefile.gnu
 GNU + CUDA 環境向けのMakefile。
-- hacapk.h
 H行列の定義。
-- loadmatrix.cpp
 hdumpで出力したバイナリファイルを読み込む機能。
-- hmvm.cpp
 main関数。
-- hmvm_seq.cpp, hmvm_seq.h
 逐次HMVMコード。
-- hmvm_omp.cpp, hmvm_omp.h
 OpenMP版HMVMコード。
-- hmvm_cuda.cu, hmvm_cuda.h
 CUDA版HMVMコード。
-- hmvm_cuda_kernels.cu, hmvm_cuda_kernels.h
 CUDAカーネル。


