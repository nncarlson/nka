# Makefile for building the test program nka_f95_test.
# 
# Modify the following lines to suit your Fortran compiler.
# "make DEBUG=1" will build a debug version; otherwise it's optimized.

# Your Fortran compiler and flags to use in optimize and debug modes.
F90 = nagfor
OPT_FLAGS = -O3 -Ounsafe -f2003 #-mtrace
DBG_FLAGS = -g -gline -f2003

# You shouldn't need to change anything below #################################

ifdef DEBUG
  F90FLAGS = $(DBG_FLAGS) $(DEFINES)
else
  DEFINES += -DNDEBUG
  F90FLAGS = $(OPT_FLAGS) $(DEFINES)
endif

PROGS= nka_f95_test
OBJ= f90_assert.o nka_type.o nka_f95_test.o

nka_f95_test: $(OBJ)
	$(F90) $(F90FLAGS) -o $@ $(OBJ)

clean:
	rm -rf *.mod *.MOD *.o
cleanest: clean
	rm $(PROGS)

.SUFFIXES: .f90 .F90
.F90.o:
	$(F90) $(F90FLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<
.f90:
	$(F90) $(F90FLAGS) $(LDFLAGS) $< -o $@

%.o : %.mod
