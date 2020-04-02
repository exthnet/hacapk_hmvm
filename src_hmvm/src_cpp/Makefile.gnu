# -*- makefile -*-

CPP  = g++ -std=c++11
CCFLAGS  += -O3 -fopenmp


######################
# Object files
OBJS  = loadmatrix.o hmvm.o
HOBJS = hacapk_cpp.hpp hmvm_seq.hpp hmvm_omp.hpp

######################
# Compile cmmands
.SUFFIXES:
.SUFFIXES: .o .cpp

$(TARGET): $(OBJS)
	$(LINK) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

.cpp.o: *.cpp
	$(CPP) -c $(CCFLAGS)  $(INCS) $<
.f90.o: *.f90
	$(F90) -c $(F90FLAGS) $(INCS) $<

all: hmvm_cpu

hmvm_cpu: $(OBJS) $(HOBJS)
	$(CPP) -o hmvm_cpu $(OBJS) -lm $(CCFLAGS)

clean:
	rm -f *.o *.mod *~
distclean:
	rm -f *.o *.mod $(TARGET) *~

rmod:
	rm -f m_*.o *.mod *~

