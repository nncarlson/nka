DEBUG=1
include ../config/config
VPATH=.:../src

OBJ= f90_assert.o fpa_kinds.o fixed_point_accelerator.o fpa_unit_test.o
unit_test: $(OBJ)
	$(F90) $(F90FLAGS) -o $@ $(OBJ)

clean:
	rm -rf $(CLEANFILES) unit_test
