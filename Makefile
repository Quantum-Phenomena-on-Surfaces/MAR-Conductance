FFLAGS =  -O  #-fopenmp 
FF = -c
FC = gfortran

LIBBLAS = -lblas 
LIBLAPACK = -llapack

LIBS = $(LIBFFTW3) $(LIBLAPACK) $(LIBBLAS)

OBJ= declarations.o algebra.o SelfEnergies.o GreenFunctions.o current.o Tools.o MAR_02.o

MAR_02.out: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o MAR_02.out

declarations.o: declarations.f90
	$(FC) $(FFLAGS) $(FF)  declarations.f90

algebra.o: algebra.f90
	$(FC) $(FFLAGS) $(FF)  algebra.f90

SelfEnergies.o: SelfEnergies.f90
	$(FC) $(FFLAGS) $(FF) SelfEnergies.f90

GreenFunctions.o: GreenFunctions.f90
	$(FC) $(FFLAGS) $(FF) GreenFunctions.f90

current.o: current.f90
	$(FC) $(FFLAGS) $(FF) current.f90

Tools.o: Tools.f90
	$(FC) $(FFLAGS) $(FF) Tools.f90

MAR_02.o: MAR_02.f90
	$(FC) $(FFLAGS) $(FF) MAR_02.f90

clean:
	@echo "CLEANING up!"
	rm -f *.o *.mod

