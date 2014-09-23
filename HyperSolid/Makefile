.SUFFIXES: .f90
LIBS = /usr/lib/liblapack.dylib \
	/usr/lib/libblas.dylib
FUN =  r_common.f90 run_variables.f90 meshgen_solid.f90 solid_variables.f90  \
	mesh_convert_variables.f90\
	s_gdof.f90 svshp.f90 getB.f90 defgrad.f90 mooney.f90 s_pen.f90 \
	s_int.f90 s_dam.f90 s_kin.f90 getM.f90 s_ess.f90 convert_abaqus_inp.f90\
	paraout.f90 r_element.f90 outnormal_2d.f90 apply_2ndbc_solid2d.f90 solve_solid.f90 main.f90

OBJ = r_common.o run_variables.o meshgen_solid.o solid_variables.o  \
	mesh_convert_variables.o\
	s_gdof.o svshp.o getB.o defgrad.o mooney.o s_pen.o \
	s_int.o s_dam.o s_kin.o getM.o s_ess.o convert_abaqus_inp.o \
	paraout.o r_element.o outnormal_2d.o apply_2ndbc_solid2d.o solve_solid.o main.o
main: .o
	gfortran -o main $(OBJ) $(LIBS)
.o:
	gfortran -c $(FUN)
clean:
	rm -rf *.vtk *.out *.o *.in *.mod main
