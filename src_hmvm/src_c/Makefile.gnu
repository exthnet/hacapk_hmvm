# -*- makefile -*-

CC   = gcc
CCFLAGS  += -Wall -O3 -fopenmp


######################
# Object files
OBJS  = loadmatrix.o hmvm.o hmvm_seq.o hmvm_omp.o


######################
# Compile cmmands
.SUFFIXES: .o .c

$(TARGET): $(OBJS)
	$(LINK) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

.c.o: *.c
	$(CC)  -c $(CCFLAGS)  $(INCS) $<
.f90.o: *.f90
	$(F90) -c $(F90FLAGS) $(INCS) $<

all: hmvm_cpu

hmvm_cpu: $(OBJS)
	$(CC) -o hmvm_cpu $(OBJS) -lm $(CCFLAGS)

clean:
	rm -f *.o *.mod *~
distclean:
	rm -f *.o *.mod $(TARGET) *~

rmod:
	rm -f m_*.o *.mod *~

