#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>



void blacs_get_(int*, int*, int*);
void blacs_pinfo_(int*, int*);
void blacs_gridinit_(int*, char*, int*, int*);
void blacs_gridinfo_(int*, int*, int*, int*, int*);
void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);

void pdgesp_( int*, int*, double*, int*, int*, int*, int*, double*, int*, int*, int*, int*);

void blacs_gridexit_(int*);
int numroc_(int*, int*, int*, int*, int*);

int main(int argc, char** argv){
	int process_rank, size_of_cluster;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);


	int n = 15;       // (Global) Matrix size

	int nprow = 1;   // Number of proc - rows
	int npcol = 1;   // Number of proc - cols
//	char uplo='L';   // Matrix is lower triangular
	char layout='R'; // Block cyclic, Row major processor mapping



	printf("my rank is %d\n", process_rank);

/*	die Matrix, die invertiert werden soll ist -I:
 *	Hierdurch sieht das Gleichungssystem wie folgt aus: 
 *	A x = b; und die Lösung ist x = -b;
 */

	double *A = malloc( n*n *sizeof(double));
	double *B = malloc( n *sizeof(double));
	double *X = malloc( n *sizeof(double));
	int *ipiv = malloc( n * sizeof(int));

	int i, j;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++)
			A[i*n + j] = i!=j ? 0 : -1; 
		B[i] = i;
		X[i] = -i;
		ipiv[i] = i;
	}

	if (process_rank == 0){
//		printf("This is process 0\n");
		printf("Size of int: %d\n", sizeof(int));
		printf("allocation and initialization done!\n");
	}   

	int zero = 0;
	int one = 1;
	int ictxt, myrow, mycol;
	int iam, nprocs;
	blacs_pinfo_(&iam, &nprocs); // BLACS rank and world size
	blacs_get_(&zero, &zero, &ictxt ); // -> Create context
	blacs_gridinit_(&ictxt, &layout, &nprow, &npcol ); // Context -> Initialize the grid
 	blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    // Compute the size of the local matrices
	int bsxA = n; //blocksize in X direction
	int bsyA = n; //blocksize in Y direction
	
	int numcA = numroc_( &n, &bsxA, &mycol, &zero, &npcol ); // number of columns stored in each process
	int numrA = numroc_( &n, &bsyA, &myrow, &zero, &nprow ); // number of rows stored in each process

	//create a descriptor for A
	int descA[9];
	int info = 0;
	int lddA = numcA > 1? numcA : 1;

	MPI_Barrier(MPI_COMM_WORLD);
	if (process_rank == 0) printf("I am process %d, \tnumrA = %d, numcA = %d, bsyA = %d, bsxA = %d\n", process_rank, numrA, numcA, bsyA, bsxA);
	MPI_Barrier(MPI_COMM_WORLD);

	descinit_( descA,  &n, &n, &bsyA, &bsxA, &zero, &zero, &ictxt, &lddA, &info);
	//descinit_( descA,  &n, &n, &one, &one, &zero, &zero, &ictxt, &lddA, &info);
	if(info != 0) {
		printf("Error in descinit A, info = %d\n", info);
		exit(1);
	}

	if (process_rank == 0) printf("Descriptor A done!\n");
	MPI_Barrier(MPI_COMM_WORLD);



	//create a desriptor for B
	
    // Compute the size of the local matrices
	int bsxB = 1; //blocksize in X direction
	int bsyB = n; //blocksize in Y direction
	
	int numcB = numroc_( &one, &bsxB, &mycol, &zero, &npcol ); // number of columns stored in each process
	int numrB = numroc_( &n, &bsyB, &myrow, &zero, &nprow ); // number of rows stored in each process
	//numcB = 1;
	int descB[9];
	int lddB = numrB > 1? numrB : 1; 	//number of local rows in b;
	


	
	MPI_Barrier(MPI_COMM_WORLD);
	if (process_rank == 0) printf("I am process %d, \tnumrB = %d, numcB = %d, bsyB = %d, bsxB = %d, \tlddB = %d\n", process_rank, numrB, numcB, bsyB, bsxB, lddB);
	MPI_Barrier(MPI_COMM_WORLD);
	
//DESCB =            1           0           4           1           4           1           0           0           4
//					   use bsxB
	descinit_( descB, &n, &one, &bsyB, &bsxB, &zero, &zero, &ictxt, &lddB, &info);
	if(info != 0) {
		printf("Error in descinit B, info = %d\n", info);
		exit(1);
	}

	if ( process_rank == 0) { 
		printf("DESCB = [");
		for (i = 0; i < 9; i++) printf(" %d", descB[i]);
		printf("]\n");
	}
	if (process_rank == 0) printf("Descriptor B done!\n");
	MPI_Barrier(MPI_COMM_WORLD);
	

	//timing
	double t0;



	t0 = MPI_Wtime();

	int nrhs = numcB;	//1;
	int ia = 0 +1;
	int ja = 0 +1; // = 0;
	int ib = 0 +1;
	int jb = 0 +1;	// = 0;
	
	if (process_rank == 0) printf("Starting PDGESV....\n");
	MPI_Barrier(MPI_COMM_WORLD);
	if (process_rank == 0) printf("I am process %d, \tia = %d, ja = %d,\tiB = %d, jB = %d\n", process_rank, ia, ja, ib, jb);
	MPI_Barrier(MPI_COMM_WORLD);
	

	//PDGESV( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
	//if (process_rank == 0) printf("n = %d, ipiv: [%d, %d, %d, %d]\n", n, ipiv[0], ipiv[1], ipiv[2], ipiv[3]);
	pdgesp_( &n, &one, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, &info);
//	pdgesv_( &n, &one, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, &info);
//	pdgesv_( 4, 1, A, ia, ja, descA, ipiv, B, ib, jb, descB, info);
	if (info!=0 && process_rank == 0) printf("ERROR DURING PDGESV, info = %d\n", info);

	if (process_rank == 0) printf("PDGESV done!\nn = %d\n", n);
	MPI_Barrier(MPI_COMM_WORLD);


	if (process_rank == 0){
	       	printf("Differenz (X - B):\n");
		for (i = 0; i < n; i++){
			printf("%f\n", X[i] - B[i]);
		}
	}
	free(A);
	free(B);
	free(X);
	MPI_Finalize();
	return 0;
}
