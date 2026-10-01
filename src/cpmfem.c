#include "functions.h"
#include <omp.h>
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//global variables (unique to MPI rank but accessible by all functions)
int myRank, size;
int my_ymin, my_ymax; // ymin and ymax rows in global coordinates
int my_nvy, my_nny, my_ndof;
int my_nvy_owned, my_nny_owned, my_ndof_owned;
int topNbr, botNbr;

// MPI comm timers
double time_cell_comm;
double time_fem_comm;

int main(int argc, char *argv[])
{
	int nrrdof, nrrdof_owned, d;
	int *dofpos;
	double *f, *u;
	double **klocal;
	int *kcol;
	double *kval;
	NOD *pn;
	VOX *pv;
	int NRc, NRc_loc, c;
	int *loc_cells; // list of local cells (stores global cell IDs)
	int v, vx, vy;
	int *csize;
	int incr, startincr;

	int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	// set num threads as env variable with export OMP_NUM_THREADS=n
	// Print out number of omp threads that the program sees
	int numthreads;
	#pragma omp parallel
	{
		numthreads = omp_get_num_threads();
		#pragma omp single
		printf("Number of omp threads on rank %d set to %d\n",myRank,numthreads);
	}

	/// INITIALIZE ///
	double start_setup = MPI_Wtime();

	// Domain Decomposition: split into y-dimension chunks
    int base = NVY / size; // we want each rank responsible for a contiguous chunk of size NVY/np
    int remainder = NVY%size; // get remainder to handle case where NVY is not evenly divisible by np
    if (myRank < size - remainder){
        my_ymin = myRank * (base);
		my_nvy_owned = base;
    } else {
        // these ranks from size - remainder to end get + 1
        my_ymin = (size - remainder) * base + (myRank - (size - remainder)) * (base + 1);
        my_nvy_owned =  base+1;
    }
    my_ymax = my_ymin + my_nvy_owned; //save min and max owned row in global coordinates 
	my_nny_owned = my_nvy_owned+1;
	
	// now define my_nvy and my_nny to include halos (2 above, 2 below)
	my_nvy = my_nvy_owned+4;
	my_nny = my_nny_owned+4;
	printf("Rank %d: global rows [%d, %d), my_nvy=%d, my_nny=%d\n",myRank, my_ymin, my_ymax, my_nvy, my_nny);

	// define top and bottom neighbors (note that first and last rank have MPI_PROC_NULL as neighbor)
	topNbr = (myRank == 0) ? MPI_PROC_NULL : myRank - 1;
	botNbr = (myRank == size - 1) ? MPI_PROC_NULL : myRank + 1;

	// Setup
   	srand(SEED); // sets global seed for pseudo-random rand()
	mt_init(); //mersenne twister using seeded rand()
   	pv = init_voxels();
	pn = init_nodes();

	// initialize cells
	NRc = init_cells(pv, &NRc_loc, &loc_cells);
	printf("Rank %d: sees %d total cells in the simulation, with %d local cells (including halo) \n",myRank,NRc, NRc_loc);
	write_cells(pv,0);

	// initialize volume (num of voxels) of each cell
	// NOTE: at initial placement each cell is only one voxel large, so no cell is spanning task boundaries yet
	csize = calloc(NRc, sizeof(int)); // keep cell volumes for all global cells (cheap memory, removes need for MPI comm)
	for(vy=2; vy<my_nvy-2; vy++) //owned voxels
   	for(vx=0; vx<NVX; vx++)
   	{
   		v = vx + vy*NVX;
		if(pv[v].ctag) { //if voxel v belongs to a cell (using global cell index), add it to that cell's volume
			csize[pv[v].ctag-1]++;
		}
	}
	// make sure all ranks have same view of csize
	MPI_Allreduce(MPI_IN_PLACE, csize, NRc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	//initialize force to zero and set BCs
	set_forces(pn);
	set_restrictions(pn);

	// local K matrix (local to a single node)
	klocal = set_klocal();

	// K matrix for all nodes on rank (including halo?):
	my_ndof = 2*(NNX*my_nny);
	kcol = calloc(10*my_ndof,sizeof(int));
	kval = calloc(10.0*my_ndof,sizeof(double));
	assembly(kcol,kval,klocal,pv);
	dofpos = calloc(my_ndof,sizeof(int));
	nrrdof = arrange_dofpos(dofpos,pn,&nrrdof_owned);

	//print for debugging
	int tot_nrrdof_owned;
	MPI_Reduce(&nrrdof_owned,&tot_nrrdof_owned,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	if (myRank == 0){ 	printf("\nTotal nrrdof owned = %d",tot_nrrdof_owned);}

	reduce_K(kcol,kval,dofpos);
	CSRMatrix Kcsr;
	Kcsr = convert_to_csr(kcol, kval); 	//structure change for fem omp thread paralellization

	// allocate memory for f and u (original code had this inside simloop)
	f=calloc(my_ndof,sizeof(double));
	u=calloc(my_ndof,sizeof(double));

	// initialize timers
	double temp;
	double time_output = 0.0;
	double time_tracforce = 0.0; //computing traction force based on distances between pts
	double time_fem = 0.0;	
	double time_pcg = 0.0;
	double time_cpm = 0.0;
	time_cell_comm = 0.0;
	time_fem_comm = 0.0;
	
	double time_setup = MPI_Wtime() - start_setup;

	/// START SIMULATION ///
	double start_simloop =  MPI_Wtime();
	for(incr=startincr; incr<NRINC; incr++)
	{
		if (myRank == 0) printf("\nSTART INCREMENT %d",incr); fflush(stdout);

		temp = MPI_Wtime();
		write_cells(pv,incr);
		time_output += MPI_Wtime() - temp;
		//printf("\n Rank %d, done write_cells", myRank); fflush(stdout);

		temp = MPI_Wtime();
		cell_forces(pv,pn,numthreads,NRc,NRc_loc,loc_cells);
		time_tracforce += MPI_Wtime() - temp;

		// FEA part // parts of this can go out the loop depending on what changes
		temp = MPI_Wtime();
		memset(f, 0, my_ndof*sizeof(double)); //set to zero instead of calling calloc each iter
		memset(u, 0, my_ndof*sizeof(double));
		place_node_forces_in_f(pn,f);
		set_disp_of_prev_incr(pn,u);		
		double temp2 = MPI_Wtime();
		//solvePCG(kcol,kval,u,f,nrrdof);
		solvePCG_csr(&Kcsr, u, f, pn);
		//printf("\nRank %d, done solvePCG_csr", myRank); fflush(stdout);
		time_pcg += MPI_Wtime() - temp2; //profile solvePCG subroutine of FEM		
		disp_to_nodes(pn,u);
		time_fem += MPI_Wtime() - temp;
		//printf("\nRank %d, done disp_to_nodes", myRank); fflush(stdout);
		
		// exchange halo buffers for pn so that strains are up to date before CPM moves
		temp = MPI_Wtime();
		exchange_node_halos(pn);
		time_fem_comm += MPI_Wtime() - temp;

		temp = MPI_Wtime();
		if(incr%5==0)
		{
			write_pstrain(pv,pn,incr);
			//printf("\nRank %d, done write_pstrain", myRank); fflush(stdout);
		}
		time_output += MPI_Wtime() - temp;

		temp = MPI_Wtime();
		CPM_moves(pv,pn,csize);
		//printf("\nRank %d, done CPM_moves", myRank); fflush(stdout);
		time_cpm +=  MPI_Wtime() - temp;

		temp = MPI_Wtime();
		update_cell_volumes(pv, csize, NRc);
		exchange_cell_halos(pv, NRc, &NRc_loc,&loc_cells);
		//printf("\n Rank %d, done exchange_cell_halos", myRank); fflush(stdout);
		time_cell_comm += MPI_Wtime() - temp;

	}
	double time_simloop = MPI_Wtime() - start_simloop;

	/// END ///
	if (myRank == 0) printf("\nSIMULATION FINISHED!"); fflush(stdout);

	// use MPI reduction to get max timing across all ranks
	double times[9] = {time_setup,time_simloop,time_output,time_tracforce,time_fem,time_pcg,time_fem_comm, time_cpm, time_cell_comm};
	if (myRank == 0) {
		MPI_Reduce(MPI_IN_PLACE, times, 9, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	} else {
		MPI_Reduce(times, NULL, 9, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	}

	if (myRank == 0){
		// Output all timing info
		printf("\nTIMING RESULTS\n");
		printf("%-20s | %s\n", "Process", "Time (sec)");
		printf("---------------------|------\n");
		// Print the table rows using left-alignment and fixed width
		printf("%-20s | %f\n", "Setup", times[0]);
		printf("%-20s | %f\n", "Total SimLoop", times[1]);
		printf("%-20s | %f\n", "Output Writing", times[2]);
		printf("%-20s | %f\n", "Traction Force", times[3]);
		printf("%-20s | %f\n", "FEM", times[4]);
		printf("%-20s | %f\n", "PCG (part of FEM)", times[5]);
		printf("%-20s | %f\n", "MPI Node Comm (part of FEM)", times[6]);
		printf("%-20s | %f\n", "CPM", times[7]);
		printf("%-20s | %f\n", "MPI Cell Comm", times[8]);
		// write the profiling data to a .csv file in the ouput folder
		write_timing(times[0],times[1],times[2],times[3],times[4],times[5],times[6],times[7],times[8],incr);
	}

	MPI_Finalize();

	free(pv); free(pn); free(klocal); free(kcol); free(kval); free(dofpos); free(u); free(f);
	free(loc_cells);
	free(Kcsr.rowptr);
	free(Kcsr.colind);
	free(Kcsr.vals);
	free(Kcsr.diag);
	return 0;
}

