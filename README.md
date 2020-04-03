# hacapk_hmvm

## src_hdump
### 概要
HACApK 1.0.0を改造し、H行列をバイナリファイルに出力するようにしたもの。
m_HACApK_hdump.f90以外のコードの多くはHACApK 1.0.0のコードを流用している。
### 使い方
makeするとbem-bb-SCM.outファイルが生成される。
MPI + C 環境であればどこでも使える、と期待される。

./bem-bb-SCM.out inputmatrix
として実行すると、hmatrix_0.binというバイナリファイルが出力される。
出力されたバイナリファイルはsrc_hmvmの入力ファイルとして使える。

元のHACApKの都合でMPIが必要だが、このバイナリ出力プログラムが通信を行うわけではないため、
mpirun等を噛ませる必要はない。

### 構成

## src_hmvm
### 概要
src_hdumpで出力したバイナリファイルを用いてH行列ベクトル積を行うテストコード群。
### 構成と使い方
#### src_c ディレクトリ：シンプルなC言語コード、CPU(OpenMP含む)向け、倍精度。

#### src_cpp ディレクトリ
C++コード、CPU(OpenMP)とGPU(CUDA)、単精度と倍精度それぞれで計算。
Makefile.gnuを使ってmake(make -f Makefile.gnu)するとhmvm_gpuまたはhmvm_cpuが生成される。
src_hdumpのbem-bb-SCM.outで出力したhmatrix_0.binを使って
./hmvm_gpu ./hmatrix_0.bin
として実行すると計算が行われる。
引数をもう1つ追加し
./hmvm_gpu ./hmatrix_0.bin 1
とすると、各実装による計算結果がそれぞれファイルに出力される。

- Makefile.gnu
  - GNU + CUDA 環境向けのMakefile。
    先頭のGPU=1を有効化するとGPU実装が有効化される。
- hacapk.h
  - H行列の定義。
- loadmatrix.cpp
  - hdumpで出力したバイナリファイルを読み込む機能。
- hmvm.cpp
  - main関数。
- hmvm_seq.cpp, hmvm_seq.h
  - 逐次HMVMコード。
- hmvm_omp.cpp, hmvm_omp.h
  - OpenMP版HMVMコード。
- hmvm_cuda.cu, hmvm_cuda.h
  - CUDA版HMVMコード。
- hmvm_cuda_kernels.cu, hmvm_cuda_kernels.h
  - CUDAカーネル。


