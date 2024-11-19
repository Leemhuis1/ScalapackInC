#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
//#include <mkl_types.h>

//#include <scalapack.h>
//#include <blacs.h>

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


	int n = 4;       // (Global) Matrix size

	int nprow = 4;   // Number of proc - rows
	int npcol = 1;   // Number of proc - cols
//	char uplo='L';   // Matrix is lower triangular
	char layout='R'; // Block cyclic, Row major processor mapping



	printf("my rank is %d\n", process_rank);

/*	die Matrix, die invertiert werden soll sieht wie folgt aus:
 *		10105	75	146	76
 *	A = 	75	10108	132	90
 *		146	132	10266	138
 *		76	90	138	10086
 *
 *	mit einer rechten Seite:
 *		5
 *	b =	10
 *		8
 *		10
 *
 *	ergibt sich für Ax = b die Lösung:
 *
 *		0.4695
 *	~x = 	0.9674		* 1e-3
 *		0.7471
 *		0.9691
*/

	double *A = malloc( 4 *sizeof(double));
	double *B = malloc( 1 *sizeof(double));
	double *X = malloc( 1 *sizeof(double));

	if (process_rank == 0){
//		printf("This is process 0 from inside if - clause\n");
		printf("Size of int: %d\n", sizeof(int));
		A[0] = 10105;
		A[1] = 75;
		A[2] = 146;
		A[3] = 76;

		B[0] = 5;
		X[0] = 0.4695e-3;
	} else if (process_rank == 1){
//		printf("This is process 1 from inside if - clause\n");
		A[0] = 75;
		A[1] = 10108;
		A[2] = 132;
		A[3] = 90;

		B[0] = 10;
		X[0] = 0.9674e-3;
	} else if (process_rank == 2){
		A[0] = 146;
		A[1] = 132;
		A[2] = 10266;
		A[3] = 139;

		B[0] = 8;
		X[0] = 0.7471e-3;
	} else if (process_rank == 3){
		A[0] = 76;
		A[1] = 90;
		A[2] = 138;
		A[3] = 10086;

		B[0] = 10;
		X[0] = 0.9691e-3;
	}    

	MPI_Barrier(MPI_COMM_WORLD);
	if (process_rank == 0) printf("A = \t%7.1f\t%7.1f\t%7.1f\t%7.1f\n", A[0], A[1], A[2], A[0]);
	int i;
	for (i = 1; i < 4; i++){
		if (process_rank == i) printf("\t%7.1f\t%7.1f\t%7.1f\t%7.1f\n", A[0], A[1], A[2], A[0]);
		MPI_Barrier(MPI_COMM_WORLD);
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
	int bsx = 4; //blocksize in X direction
	int bsy = 1; //blocksize in Y direction
	
	int numc = numroc_( &n, &bsx, &mycol, &zero, &npcol ); // number of columns stored in each process
	int numr = numroc_( &n, &bsy, &myrow, &zero, &nprow ); // number of rows stored in each process

	MPI_Barrier(MPI_COMM_WORLD);
	for (i = 0; i < 4; i++){
		if (process_rank == i) printf("I am process %d, \tnumr = %d, numc = %d,\tcoords %d %d\n", process_rank, numr, numc, myrow, mycol);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	//create a descriptor for A
	int descA[9];
	int info = 0;
	int lddA = numc > 1? numc : 1;

	int bfr = n / numc; //blocking factor columns 	// = 1 ( = bsy)
	int bfc = n / numr; //blocking factor rows 	// = 4 ( = bsx)
		
	MPI_Barrier(MPI_COMM_WORLD);
	for (i = 0; i < 4; i++){
		if (process_rank == i) printf("I am process %d, \tnumr = %d, numc = %d,\tbfr = %d, bfc = %d\n", process_rank, numr, numc, bfr, bfc);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	descinit_( descA,  &n, &n, &bfr, &bfc, &zero, &zero, &ictxt, &lddA, &info);
	if(info != 0) {
		printf("Error in descinit A, info = %d\n", info);
		exit(1);
	}

	if (process_rank == 0) printf("Descriptor A done!\n");
	MPI_Barrier(MPI_COMM_WORLD);



	//create a desriptor for B
	int descB[9];
	int numcB = 1;
	int numrB = 1;
	int lddB = numcB > 1? numcB : 1; 	//number of local rows in b;
	

	int bfcB = 1;	//blocking factor - col (blocksize in cols)
	int bfrB = 1;	//blocking factor - row (blocksize in rows)
	
	MPI_Barrier(MPI_COMM_WORLD);
	for (i = 0; i < 4; i++){
		if (process_rank == i) printf("I am process %d, \tnumrB = %d, numcB = %d,\tbfrB = %d, bfcB = %d\n", process_rank, numrB, numcB, bfrB, bfcB);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	descinit_( descB, &n, &one, &bfrB, &bfcB, &zero, &zero, &ictxt, &lddB, &info);
	if(info != 0) {
		printf("Error in descinit B, info = %d\n", info);
		exit(1);
	}

	if (process_rank == 0) printf("Descriptor B done!\n");
	MPI_Barrier(MPI_COMM_WORLD);
	


	//timing
	double t0;



	t0 = MPI_Wtime();

	int nrhs = numcB;	//1;
	int ia = myrow;
	int ib = myrow;
	int ja = mycol; // = 0;
	int jb = mycol;	// = 0;
	int ipiv[4];
	for (i = 0; i<4; i++) ipiv[i] = i;

	if (process_rank == 0) printf("IPIV done!\n");
	if (process_rank == 0) printf("Starting PDGESV....\n");
	MPI_Barrier(MPI_COMM_WORLD);


//	n = 2;
	//PDGESV( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
	if (process_rank == 0) printf("n = %d, ipiv: [%d, %d, %d, %d]\n", n, ipiv[0], ipiv[1], ipiv[2], ipiv[3]);
	pdgesp_( &n, &one, A, &ia, &ja, descA, &zero, B, &ib, &jb, descB, &info);
//	pdgesv_( &n, &one, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, &info);
//	pdgesv_( 4, 1, A, ia, ja, descA, ipiv, B, ib, jb, descB, info);
	if (info!=0 && process_rank == 0) printf("ERROR DURING PDGESV, info = %d\n", info);

	if (process_rank == 0) printf("PDGESV done!\n");
	MPI_Barrier(MPI_COMM_WORLD);


	free(A);
	free(B);
	free(X);
	MPI_Finalize();
	return 0;
}