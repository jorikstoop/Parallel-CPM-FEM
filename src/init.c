// file: init.c
#include "functions.h"
#include <mpi.h>
#include <assert.h>
#include <stdint.h>

////////////////////////////////////////////////////////////////////////////////
VOX* init_voxels(void)
{
	VOX* pv;
	int v, vx, vy;
	int i;

	//VOX is a struct data type containing int ctag,
	//which is either the id of the occupying cell or 0 if no cell
	// use calloc to initialize all voxels as unoccupied
   	pv = calloc( my_nvy * NVX, sizeof(VOX)); 

	return pv;
}

////////////////////////////////////////////////////////////////////////////////
NOD* init_nodes(void)
{
	NOD* pn;
	int n, nx, ny;

	// NOD is a struct data type containing force, disp, and bc info
	// use calloc to set initial force and disp to zero, no BCs defined yet
   	pn = calloc( my_nny * NNX, sizeof(NOD));

	return pn;
}

////////////////////////////////////////////////////////////////////////////////
// simple hash of global voxel index + seed
uint32_t hash_vox(int gvx, int gvy, int seed) {
    uint32_t h = (uint32_t)(gvx * 2654435761u ^ gvy * 2246822519u ^ seed);
    h ^= h >> 16; h *= 0x45d9f3b; h ^= h >> 16;
    return h;
}

int init_cells(VOX* pv, int *NRC_loc, int **loc_cells)
{
	int v, vx, vy;
	int gy;
	int loc_NRc, NRc;
	int offset, loc_idx;
	double r01;
	double d; int dx, dy; // distance to center

	// 1. place cells locally with temporary local cell ID
	loc_NRc = 0; // local number of cells
	for(vy=2; vy<my_nvy-2; vy++) // loop through owned voxels
   	for(vx=0; vx<NVX; vx++)
   	{
   		v = vx + vy*NVX; //local v
		gy = my_ymin + (vy-2); // global y-position

		if((vx>0)&&(vx<NVX-1)) // exclude outer edges in x-dir
		{
			if (myRank == 0 && vy==2) // for first rank exclude top (note that local vy==2 is global y==0)
				continue;

			if (myRank == size-1 && vy==my_nvy-3) // for last rank exclude bottom (note that local vy==my_nvy-3 is global y==NNV-1)
				continue;

			//r01 = rand()/(double)RAND_MAX;
			r01 = (hash_vox(vx, gy, SEED) & 0xFFFFFF) / (double)0x1000000; // task invariant rng based on global position, useful for MPI
		
			//if (r01 < .25 / TARGETVOLUME) // place cells with 25% of max density
			// alternate placement methods:
			if((vx==NVX/2)&&(gy==NVY/2)) // to place cell at exact global center
			//if(((vx==NVX/2-7)||(vx==NVX/2+7))&&(gy==NVY/2)) // seed two cells at offset from global center
			//dx=vx-NVX/2; dy=gy-NVY/2; d=sqrt(dx*dx+dy*dy); if((d<NVX/8.0) && (r01<1.5/TARGETVOLUME)) // seed inside circular region
			{
				loc_NRc++;
				pv[v].ctag = loc_NRc;
               
			}
		}

	}

	// 2. allgather the local counts so that each rank can figure out the global cell IDs
	int *counts = malloc(size * sizeof(int));
	MPI_Allgather(&loc_NRc, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

	NRc = 0; // global cell count
	offset = 0; // start idx for this rank
    for (int r = 0; r < size; r++){
		if (r < myRank) offset += counts[r];
        NRc += counts[r];
	}
    free(counts);

	// 3. re-write ctags with global cell ID (since no cells die or divide, global ID should be static throughout simulation)
    for (vy = 2; vy < my_nvy - 2; vy++)	// loop through owned voxels
    for (vx = 0; vx < NVX; vx++)
    {
        v = vx + vy * NVX;
        if (pv[v].ctag > 0) {
            pv[v].ctag += offset; //convert to global cell ID
        }
    }

	// 4. Exchange halo information so all tasks are aware of cells near boundaries

	// rank sends first two owned rows to topNbr, and botNbr receives into its bottom two rows (the halo)
	MPI_Sendrecv(&pv[2 * NVX], 2 * NVX, MPI_INT, topNbr, 0, &pv[(my_nvy-2) * NVX], 2 * NVX, MPI_INT, botNbr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// rank sends last two owned rows to botNbr, and topNbr receives into its top two rows (the halo)
	MPI_Sendrecv(&pv[(my_nvy-4) * NVX], 2 * NVX, MPI_INT, botNbr, 1, &pv[0], 2 * NVX, MPI_INT, topNbr, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// 5. Fill local cells list after halo exchange
	// at this moment each cell is only 1 pixel, so don't need to check for uniqueness when filling loc_cells
	*loc_cells = malloc(NRc * sizeof(int)); //over allocate to avoid dynamic re-allocation
	loc_idx = 0;
    for (vy = 0; vy < my_nvy; vy++)	// loop through all local voxels
    for (vx = 0; vx < NVX; vx++)
    {
        v = vx + vy * NVX;
        if (pv[v].ctag > 0) {
            (*loc_cells)[loc_idx] = pv[v].ctag; //store global cell ID in local list
			loc_idx++;
        }
    }
	*NRC_loc = loc_idx; // save local cells after halo exchange

	return NRc;
}





////////////////////////////////////////////////////////////////////////////////
void set_forces(NOD* pn)
{
   	int n, nx, ny;
	int global_ny;

	double a = (0.0/6.0) * 3.1416;

	// reset all local nodal forces to zero (including halo)
	for(n=0; n< NNX * my_nny; n++)
	{
		pn[n].fx = .0;
		pn[n].fy = .0;
	}

	// set boundary conditions
   	for(ny=2; ny<my_nny-2; ny++)
   	for(nx=0; nx<NNX; nx++)
   	{
   		n = nx + ny*NNX;
		global_ny = my_ymin + (ny - 2);

		// top plate (iy==0) loading (only)
		if(global_ny == 0)
		{
			assert(myRank==0); //sanity check, should only happen on first rank
        	pn[n].fx +=  sin(a)*cos(a)*FORCE;
			pn[n].fy += -cos(a)*cos(a)*FORCE;
       	}
		// bot plate (iy==NNY-1) loading
		if(global_ny == NNY-1)
       	{
			assert(myRank == size-1);
        	pn[n].fx += -sin(a)*cos(a)*FORCE;
			pn[n].fy +=  cos(a)*cos(a)*FORCE;
		}
		// left plate (ix==0) loading
		if(nx==0)
		{
        	pn[n].fx += -sin(a)*sin(a)*FORCE;
			pn[n].fy +=  sin(a)*cos(a)*FORCE;
		}
		// right plate (ix==NNX-1) loading
		if(nx==NNX-1)
		{
        	pn[n].fx +=  sin(a)*sin(a)*FORCE;
			pn[n].fy += -sin(a)*cos(a)*FORCE;
        }
	}

   	for(ny=2; ny<my_nny-2; ny++)
   	for(nx=0; nx<NNX; nx++)
   	{
   		n = nx + ny*NNX;
		global_ny = my_ymin + (ny - 2);

      	// for loading on the side of a plate, forces are lower
		if((nx==0)||(nx==NNX-1)||(global_ny==0)||(global_ny==NNY-1))
		{
			pn[n].fx *= .5;
         	pn[n].fy *= .5;
		}
   	}
}

////////////////////////////////////////////////////////////////////////////////
void set_restrictions(NOD* pn)
{
	int n, nx, ny;
	int global_ny;
	
	//mark all boundary nodes
	for(ny=2; ny<my_nny-2; ny++)
   	for(nx=0; nx<NNX; nx++)
   	{
   		n = nx + ny*NNX;
		global_ny = my_ymin + (ny - 2);

		// if at x boundary
		// or at top for first rank
		// or at bottom for last rank
		if((nx==0)||(nx==NNX-1)||(global_ny==0)||(global_ny==NNY-1))
		{
			pn[n].restrictx=TRUE;
			pn[n].restricty=TRUE;
		}
	}

}

int is_owned_node(int global_ny)
{
    if (myRank == size - 1)
        return (global_ny >= my_ymin && global_ny <= my_ymax);
    else
        return (global_ny >= my_ymin && global_ny < my_ymax);
}