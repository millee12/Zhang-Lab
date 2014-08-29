.SUFFIXES: .f90
LIBS = /usr/lib/liblapack.dylib \
	/usr/lib/libblas.dylib
OBJ = init_solid.o meshgen_solid.o solid_variables.o  \
	mesh_convert_variables.o\
	s_gdof.o svshp.o getB.o defgrad.o mooney.o s_pen.o \
	s_int.o s_kin.o getM.o s_ess.o convert_abaqus_inp.o paraout.o main.o
main: $(OBJ)
	gfortran -o main $(OBJ) $(LIBS)
.f90.o:
	gfortran -c $<
cleanall:
	rm -rf *.vtk *.out *.o *.in main
clean:
	rm -rf fem* *.out *.o *.mod main
