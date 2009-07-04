#DEBUG=1
include ../config/config
VPATH=.:../src

PROGS= unit_test
OBJ= f90_assert.o nka_kinds.o nonlinear_krylov_accelerator.o nka_unit_test.o

unit_test: $(OBJ)
	$(F90) $(F90FLAGS) -o $@ $(OBJ)

clean:
	rm -rf $(CLEANFILES) unit_test
