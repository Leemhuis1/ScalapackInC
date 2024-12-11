#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>

//#include <mkl_types.h>

//#include <scalapack.h>
//#include <blacs.h>

void blacs_get_(int*, int*, int*);
void blacs_pinfo_(int*, int*);
void blacs_gridinit_(int*, char*, int*, int*);
void blacs_gridinfo_(int*, int*, int*, int*, int*);
void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);

void pdgesv_( int*, int*, double*, int*, int*, int*, int*, double*, int*, int*, int*, int*);

void blacs_gridexit_(int*);
int numroc_(int*, int*, int*, int*, int*);

int main(int argc, char** argv){
	int process_rank, size_of_cluster;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);


	int n = 6;       // (Global) Matrix size
	int m = n;
	int nrhs = 1;

	int nprow = 2;   // Number of proc - rows
	int npcol = 2;   // Number of proc - cols
	char layout='R'; // Block cyclic, Row major processor mapping

	
	// Blas content and grid setup
	int ictxt, myrow, mycol;
	int zero = 0;
	int one = 1;
	int info = 0;

	blacs_pinfo_(&process_rank, &size_of_cluster); // BLACS rank and world size (NOT NEEDED (?))
	blacs_get_(&zero, &zero, &ictxt ); // -> Create context
	blacs_gridinit_(&ictxt, &layout, &nprow, &npcol ); // Context -> Initialize the grid
	blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    // Block Sizes
	int bsx = 2; //blocksize in X direction
    	int bsy = 2; //blocksize in Y direction
    
    // Computing local matrix sizes
    	int numc = numroc_( &n, &bsx, &mycol, &zero, &npcol ); // number of columns stored in each process
	int numr = numroc_( &m, &bsy, &myrow, &zero, &nprow ); // number of rows stored in each process

    // Leading Dims
    	int lddA = numr > 1? numr : 1;
    	int lddB = numr > 1? numr : 1; 	//number of local rows in b;
    
    
    // Local arrays
    	double *A = malloc( numc * numr *sizeof(double));
	double *B = malloc( numr *sizeof(double));
	double *X = malloc( numr *sizeof(double));

    
	//setting local A
		
	//**********************************************
	//!!!! data storage column wise for FORTRAN !!!! 
	//**********************************************

	if (myrow == 0 && mycol == 0) {
		printf("myrank: %d, numr: %d, numc:%d\n", process_rank, numr, numc);
		for (int j = 0; j < numc; j++){
			for (int i = 0; i < numr; i++){
				A[j * numr + i] = i > j? 0:1;
			}
		}	
	}
 	sleep(0.05);
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize output
	if (myrow == 0 && mycol == 1){
		printf("myrank: %d, numr: %d, numc:%d\n", process_rank, numr, numc);
		A[0] = 1; A[4] = 1;
		A[1] = 1; A[5] = 1;
		A[2] = 0; A[6] = 0;
		A[3] = 0; A[7] = 0;
	}
	sleep(0.05);
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize output
	if (myrow == 1 && mycol == 0){
		printf("myrank: %d, numr: %d, numc:%d\n", process_rank, numr, numc);
		A[0] = 0; A[2] = 0; A[4] = 1; A[6] = 1;
		A[1] = 0; A[3] = 0; A[5] = 1; A[7] = 1;
	}
	sleep(0.05);
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize output
	if (myrow == 1 && mycol == 1){
		printf("myrank: %d, numr: %d, numc:%d\n", process_rank, numr, numc);
		A[0] = 1; A[2] = 1; 
		A[1] = 0; A[3] = 1;
	}

	sleep(0.05);
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize output

/*
	if (myrow == mycol) { //"I am a process on the diagonal on the process grid"
		for (int j = 0; j < numc; j++)
			for (int i = 0; i < numr; i++)
				A[j * numr + i] = i > j? 0: 1; //if row ind. > col ind put 0 else 1
	} else if( myrow > mycol) { // "I am a process under the diagonal on the process grid"
		for (int i = 0; i < numc * numr; i++) A[i] = 0; //setting everthing to 0
	} else if( myrow < mycol) { // "I am a process above the diagonal on the process grid"
		for (int i = 0; i < numc * numr; i++) A[i] = 1; //setting everthing to 1
       	}
*/



	//setting local B //only on first process row
	if (mycol == 0 && myrow == 0){
		B[0] = 6;
		B[1] = 5;
		B[2] = 2;
		B[3] = 1;
	} else if (mycol == 0 && myrow == 1){
		B[0] = 4;
		B[1] = 3;
	}
	/*
	if (mycol == 0)
		for (int i = 0; i < numr; i++)
			B[i] = m - numr*myrow - i;
*/

	//printing B
	for (int i = 0; i < size_of_cluster; i++) {
        	if (process_rank == i && mycol == 0) { //only on first process column
            		for (int j = 0; j< numr; j++){
                		printf("rank %d \t B[%d] = %f \n", process_rank, j, B[j]);
		                fflush(stdout);
            		}
	        }
        	sleep(0.05);
	        MPI_Barrier(MPI_COMM_WORLD); // Synchronize output
	}

	//printing A
	MPI_Barrier(MPI_COMM_WORLD);
	if (process_rank == 0) printf("A = \n");

	for (int k = 0; k < size_of_cluster; k++){
		if (process_rank == k) {
			for (int j = 0; j < numr; j++) {
				for (int i = 0; i < numc; i++) printf("\t%4.1f", A[i*numr + j]);
				printf("\t held by %d \n", process_rank);
			}
			printf("\n");
		}
		MPI_Barrier(MPI_COMM_WORLD); sleep(0.05);
	}


    
	
	for (int i = 0; i < 4; i++){
		if (process_rank == i) printf("I am process %d, \tnumr = %d, numc = %d,\tcoords %d %d\n", process_rank, numr, numc, myrow, mycol);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	//create a descriptor for A and B
	int descA[9], descB[9];
	

	descinit_( descA, &m, &n, &bsy, &bsx, &zero, &zero, &ictxt, &lddA, &info);
	if(info != 0) {
		printf("Error in descinit A, info = %d\n", info);
		exit(1);
	}
	

	descinit_( descB, &m, &nrhs, &bsx, &one, &zero, &zero, &ictxt, &lddB, &info);
	if(info != 0) {
		printf("Error in descinit B, info = %d\n", info);
		exit(1);
	}
	
    // pivot array
	int ipiv[numr +bsx ]; //[4];
	int ia = 1,  ja = 1; 
    	int ib = 1,  jb = 1;	
    
	//timing
	double t0;

	t0 = MPI_Wtime();

	
	
	if (process_rank == 0) printf("Starting PDGESV....\n");
	
	if (process_rank == 0) printf("n = %d, ipiv: [%d, %d, %d, %d]\n", n, ipiv[0], ipiv[1], ipiv[2], ipiv[3]);
	pdgesv_( &n, &one, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, &info);
	if (info!=0 && process_rank == 0) printf("ERROR DURING PDGESV, info = %d\n", info);

	if (process_rank == 0) printf("PDGESV done!\n");
	MPI_Barrier(MPI_COMM_WORLD);

    	for (int i = 0; i < size_of_cluster; i = i+2) {
        	if (process_rank == i) {
            		for (int j = 0; j< numr; j++){
                		printf("rank %d \t B[%d] = %e \n", process_rank, j, B[j]);
                		fflush(stdout);
            		}
        	}
        	sleep(0.05);
        	MPI_Barrier(MPI_COMM_WORLD); // Synchronize output
	}
	
	//TODO: gather B
 


	free(A);
	free(B);
	free(X);
	MPI_Finalize();
	return 0;
}
