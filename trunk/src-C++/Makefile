# Edit this to point to your installation of Trilinos.
TRILINOS_ROOT=/home/berndt/Trilinos/trilinos-10.2.0-serial

# There is no need to edit below this line.
include $(TRILINOS_ROOT)/include/Makefile.export.NOX

CXX = $(NOX_CXX_COMPILER)
CXXFLAGS = $(NOX_CXX_FLAGS) -I$(TRILINOS_ROOT)/include


all: nlkain-trilinos.exe

NKA.o: NKA.cpp NKA.H
NKADirection.o: NKADirection.cpp NKADirection.H NKA.H
NKADirFactory.o: NKADirFactory.cpp NKADirFactory.H NKADirection.H

nlkain-trilinos.o: nlkain-trilinos.cpp NKA.H NKADirection.H NKADirFactory.H



nlkain-trilinos.exe: nlkain-trilinos.o NKA.o NKADirection.o NKADirFactory.o
	$(CXX)  NKADirFactory.o NKADirection.o NKA.o \
		nlkain-trilinos.o -o nlkain-trilinos.exe \
		-L $(TRILINOS_ROOT)/lib $(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES)

clean:
	rm -f *.o nlkain-trilinos.exe
