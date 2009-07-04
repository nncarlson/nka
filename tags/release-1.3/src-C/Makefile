# Makefile for building the test program nka_c_test.
#
# Modify the following lines to suite your C compiler
# "make DEBUG=1" builds a debug version; otherwise it's optimized.

# Your C compiler and flags to use in optimize and debug modes.
CC = gcc
OPT_FLAGS = -O3
DBG_FLAGS = -g

# You shouldn't need to change anything below #################################

ifdef DEBUG
  CFLAGS = $(DBG_FLAGS) $(DEFINES)
else
  DEFINES = -DNDEBUG
  CFLAGS = $(OPT_FLAGS) $(DEFINES)
endif

PROGS = nka_c_test
OBJ = nonlinear_krylov_accelerator.o nka_c_test.o

nka_c_test: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

clean:
	rm -rf *.o
cleanest: clean
	rm -rf $(PROGS)
