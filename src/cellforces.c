// file: cellforces.c
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h> 
#include <assert.h>
////////////////////////////////////////////////////////////////////////////////
void cell_forces(VOX* pv, NOD* pn, int nthreads, int NRc, int NRc_loc, int* loc_cells)
{
	int idx, c;
	int n,nx,ny;
	int global_ny, owned;
	int v,vx,vy, cnttag;    
	double xi, yi;
	int NN_loc = my_nny * NNX;

	// allocate local and global arrays to enable partial calculation of cell centroids by each rank
    int    *local_count  = calloc(NRc, sizeof(int));
    double *local_sumx   = calloc(NRc, sizeof(double));
    double *local_sumy   = calloc(NRc, sizeof(double));
    int    *global_count = calloc(NRc, sizeof(int));
    double *global_sumx  = calloc(NRc, sizeof(double));
    double *global_sumy  = calloc(NRc, sizeof(double));

	//1. accumulate local sums for centroid calculation
	int *count_private  = calloc(nthreads * NRc, sizeof(int)); //thread-safe arrays
	double *sumx_private = calloc(nthreads * NRc, sizeof(double));
	double *sumy_private = calloc(nthreads * NRc, sizeof(double));
	#pragma omp parallel private(idx, c, n, nx, ny, global_ny, owned, v, vx, vy, cnttag, xi, yi)
	{
		int tid = omp_get_thread_num();
		int    *count  = &count_private[tid * NRc];
		double *sumx   = &sumx_private[tid * NRc];
		double *sumy   = &sumy_private[tid * NRc];

		#pragma omp for schedule(dynamic)
		for(ny=0; ny<my_nny; ny++) //loop through local nodes
   		for(nx=1; nx<NNX-1; nx++) // exclude boundary nodes
   		{
   			n = nx + ny*NNX;
		
			global_ny = my_ymin + (ny - 2);
			//if (nx == NNX/2) printf("\nRank %d, global_ny, my_ymin, my_ymax = %d, %d, %d",myRank,global_ny, my_ymin, my_ymax);
			owned = is_owned_node(global_ny);

            // skip if not owned
            if (!owned) continue;

            // exclude boundary nodes
            if (global_ny == 0 || global_ny == NNY-1) continue;

			for(int i=0;i< NRc_loc;i++) //check all local cells
        	{
			c = loc_cells[i]; // get global cell ID
			assert(c > 0 && c <= NRc);
			idx = c-1;
			cnttag = 0;
			// check all four voxels surrounding node n
			for(vy=ny-1; vy<ny+1; vy++)
			for(vx=nx-1; vx<nx+1; vx++)
			{
				v = vx + vy*NVX;
				// check if voxel belongs to cell c
				if(pv[v].ctag == c)
					cnttag++;
			}
			if(cnttag>0) // all cell nodes
			{
				xi = nx * VOXSIZE;
				//yi = ny * VOXSIZE;
				yi = (my_ymin + (ny - 2)) * VOXSIZE;  // global y position
				count[idx]++;
				sumx[idx] += xi;
				sumy[idx] += yi;
			}
			}
		}
	}//end parallel
	//merge
	for (int t = 0; t < nthreads; t++) {
		for (idx = 0; idx < NRc; idx++) {
			local_count[idx] += count_private[t * NRc + idx];
			local_sumx[idx]  += sumx_private[t * NRc + idx];
			local_sumy[idx]  += sumy_private[t * NRc + idx];
		}
	}
    free(count_private);
    free(sumx_private);
	free(sumy_private);

	// reduce to get global counts and sums
    MPI_Allreduce(local_count, global_count, NRc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_sumx, global_sumx, NRc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_sumy, global_sumy, NRc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// //sanity check
	// int total_nodes = 0;
	// for (int c=0; c<NRc; c++)
	// 	total_nodes += global_count[c];
	// printf("\n Rank %d, Total nodes counted = %d\n", myRank, total_nodes);

	#pragma omp parallel for private(idx, c, n, nx, ny, global_ny, v, vx, vy, cnttag, xi, yi) schedule(dynamic)
    for(ny=0; ny<my_nny; ny++) //loop over local nodes
    for(nx=1; nx<NNX-1; nx++) //exlcude boundary nodes
    {
        n = nx + ny*NNX;

		global_ny = my_ymin + (ny - 2);
		owned = is_owned_node(global_ny);

        // skip if not owned (only apply force to owned nodes)
        if (!owned) continue;

        // exclude boundary nodes
        if (global_ny == 0 || global_ny == NNY-1) continue;

        pn[n].fx = 0.0;
        pn[n].fy = 0.0;

		for(int i=0;i<NRc_loc;i++) //check all local cells
        {
			c = loc_cells[i]; // get global cell ID
			assert(c > 0 && c <= NRc);
            idx = c-1;
			cnttag = 0;
            for(vy=ny-1; vy<ny+1; vy++)
            for(vx=nx-1; vx<nx+1; vx++)
            {
                v = vx + vy*NVX;
                if(pv[v].ctag == c)
                    cnttag++;
            }

            if(cnttag>0 && global_count[idx] > 0)
            {
				xi = nx * VOXSIZE;
				//yi = ny * VOXSIZE;
				yi = (my_ymin + (ny - 2)) * VOXSIZE;  // global y position
                pn[n].fx += CELLFORCE * (global_sumx[idx] - global_count[idx] * xi); // n is unique per thread, no race condition
                pn[n].fy += CELLFORCE * (global_sumy[idx] - global_count[idx] * yi);
            }
        }
    }

	// free temp arrays
    free(local_count);
    free(local_sumx);
    free(local_sumy);
    free(global_count);
    free(global_sumx);
    free(global_sumy);

}

////////////////////////////////////////////////////////////////////////////////
