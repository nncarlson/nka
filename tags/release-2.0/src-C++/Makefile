# Edit this to point to your installation of Trilinos.
TRILINOS_INSTALL=/backup/berndt/devel/tpls/trilinos-11.0.3
BOOST_INSTALL=/backup/berndt/devel/tpls

# There is no need to edit below this line.
include $(TRILINOS_INSTALL)/include/Makefile.export.NOX

CXX = $(NOX_CXX_COMPILER)
CXXFLAGS = $(NOX_CXX_FLAGS) -I$(TRILINOS_INSTALL)/include -I$(BOOST_INSTALL)/include


all: nlkain-trilinos.exe

NKA.o: NKA.cpp NKA.H
NKADirection.o: NKADirection.cpp NKADirection.H NKA.H
NKADirFactory.o: NKADirFactory.cpp NKADirFactory.H NKADirection.H

BroydenSMDirection.o: BroydenSMDirection.cpp BroydenSMDirection.hpp
BroydenSMDirFactory.o: BroydenSMDirFactory.cpp BroydenSMDirFactory.hpp BroydenSMDirection.hpp


nlkain-trilinos.o: nlkain-trilinos.cpp NKA.H NKADirection.H NKADirFactory.H BroydenSMDirection.hpp BroydenSMDirFactory.hpp



nlkain-trilinos.exe: nlkain-trilinos.o NKA.o NKADirection.o NKADirFactory.o BroydenSMDirection.o BroydenSMDirFactory.o
	$(CXX)  NKADirFactory.o NKADirection.o NKA.o BroydenSMDirection.o BroydenSMDirFactory.o\
		nlkain-trilinos.o -o nlkain-trilinos.exe \
		-L $(TRILINOS_INSTALL)/lib $(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES)

clean:
	rm -f *.o nlkain-trilinos.exe
