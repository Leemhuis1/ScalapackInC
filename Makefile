CC = mpicc
# CFLAGS = -O3  
SCALAPACKLIB = -L/home/wmaisrv1/leemhuis/Downloads/scalapack-2.2.0 
LAPACKLIB = -llapack
BLASLIB = -lblas

OBJ = large_matrix.o hello_world.o
LIBS = -lm -ldl -lgfortran $(LAPACKLIB) $(BLASLIB) -lscalapack

############################################################################
#
#  C preprocessor definitions:  set CDEFS to one of the following:
#
#     -DNoChange (fortran subprogram names are lower case without any suffix)
#     -DUpCase   (fortran subprogram names are upper case without any suffix)
#     -DAdd_     (fortran subprogram names are lower case with "_" appended)

CDEFS = -DAdd_

all: hello_world large_matrix

%.o: %.c
	$(CC) -c -o $@ $< $(CDEFS) $(CFLAGS) $(SCALAPACKLIB)

large_matrix: large_matrix.o #$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(SCALAPACKLIB) $(LIBS)

hello_world: hello_world.o #$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(SCALAPACKLIB) $(LIBS)

.PHONY: clean

clean: 
	rm -f -v *.o hello_world large_matrix
