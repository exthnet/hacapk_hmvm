# -*- makefile -*-

# comment out if you don't want to use GPU
GPU = 1

ifeq ($(GPU),1)
CPP  = g++ -std=c++11
CCFLAGS  = -O3 -fopenmp -D_USE_CUDA
NVCC = /usr/local/cuda/bin/nvcc
NVCCFLAGS = -Xcompiler "-O3 -fopenmp" \
--generate-code arch=compute_53,code=sm_53 \
--generate-code arch=compute_60,code=sm_60 \
--generate-code arch=compute_70,code=sm_70
# 53=JetsonTX1, 60=Pascal, 70=Volta
TARGET = hmvm_gpu
else
CPP  = g++ -std=c++11
CCFLAGS  = -O3 -fopenmp
TARGET = hmvm_cpu
endif

######################
# Object files
OBJS  = loadmatrix.o hmvm.o hmvm_omp.o hmvm_seq.o
CUOBJS = hmvm_cuda.o hmvm_cuda_kernels.o
HOBJS = hacapk.h

######################
# Compile cmmands
.SUFFIXES:
.SUFFIXES: .o .cpp
.SUFFIXES: .o .cu

$(TARGET): $(OBJS)
	$(LINK) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

.cpp.o: *.cpp $(HOBJS)
	$(CPP) -c $(CCFLAGS)  $(INCS) $<
.cu.o: *.cu
	$(NVCC) -c $(NVCCFLAGS) $(INCS) $<

all: $(TARGET)

hmvm_cpu: $(OBJS) $(HOBJS)
	$(CPP) -o hmvm_cpu $(OBJS) -lm $(CCFLAGS)
hmvm_gpu: $(OBJS) $(HOBJS) $(CUOBJS)
	$(NVCC) -o hmvm_gpu $(OBJS) $(CUOBJS) -lm $(NVCCFLAGS)

clean:
	rm -f *.o *.mod *~
distclean:
	rm -f *.o *.mod $(TARGET) *~

rmod:
	rm -f m_*.o *.mod *~

