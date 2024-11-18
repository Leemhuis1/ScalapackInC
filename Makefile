CC = mpicc
CFLAGS = -O3  
SCALAPACKLIB = -L/home/wmaisrv1/leemhuis/Downloads/scalapack-2.2.0 
LAPACKLIB = -llapack
BLASLIB = -lblas

OBJ = hello_world.o
LIBS = -lm -ldl -lgfortran $(LAPACKLIB) $(BLASLIB) -lscalapack

############################################################################
#
#  C preprocessor definitions:  set CDEFS to one of the following:
#
#     -DNoChange (fortran subprogram names are lower case without any suffix)
#     -DUpCase   (fortran subprogram names are upper case without any suffix)
#     -DAdd_     (fortran subprogram names are lower case with "_" appended)

CDEFS = -DAdd_



%.o: %.c
	$(CC) -c -o $@ $< $(CDEFS) $(CFLAGS) $(SCALAPACKLIB)

hello_world: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(SCALAPACKLIB) $(LIBS)

.PHONY: clean

clean: 
	rm -f -v *.o hello_world
