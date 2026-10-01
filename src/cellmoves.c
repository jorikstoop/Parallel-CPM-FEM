// file: cellmoves.c
#include "functions.h"
#include <mpi.h> 
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
void CPM_moves(VOX* pv, NOD* pn, int* csize)
// cellular potts model: one Monte Carlo step
{
	int i;
	int NRsteps = my_nvy * NVX;
	int xs, xt; // source and target pixel
	int xtx,xty; // x and y position of target pixel
	int ttag, stag; // target and source label
	int nbs[8],pick; // neighbors of target pixel
	BOOL go_on;
	double dH, prob;

	for(i=0;i<NRsteps;i++) //loop through all local voxels
	{
		//xt = (rand()*NV/RAND_MAX); // pick random element

		//NOTE: because this is different than serial RNG it will not be exactly comparable
		xt = mt_random()%NRsteps; // pick random element (LOCAL to rank)
		xty = xt/NVX;
		xtx = xt%NVX;
		int global_xty = my_ymin + (xty-2);

		if((xtx>0) && (xtx<NVX-1) && (global_xty>0) && (global_xty<NVY-1)) // exclude outer rim
		{
			//only check owned 
			if (xty < 2 || xty >= my_nvy-2) continue; //exclude halo region (only update owned voxels, which requires halo info)

			nbs[0]=xt-1+NVX; nbs[1]=xt+NVX; nbs[2]=xt+1+NVX;
			nbs[7]=xt-1;                    nbs[3]=xt+1;
			nbs[6]=xt-1-NVX; nbs[5]=xt-NVX; nbs[4]=xt+1-NVX;
			pick = mt_random()%8;
			xs = nbs[pick]; // pick random neighbor

			ttag = pv[xt].ctag; //global cell coord of target
			stag = pv[xs].ctag;

			go_on = 0;
			if(ttag!=stag) //don't bother if no difference
			{
        		go_on = 1;
        		if(ttag) // if a cell in xt (retracting)
				{
            		//if (splitcheckCCR(pv,csize,xt,ttag)) //make sure not to split cell
					if (splitcheck_local(pv,xt,ttag))
                		go_on = 0;
            		if(csize[ttag-1]==1) // cell cannot disappear (constraint may be removed)
                		go_on = 0;
				}
			}

			if(go_on)
			{
				// calculate probability based on hamiltonian function
        		dH = calcdH(pv,pn,csize,xt,xs,pick,ttag,stag);
        		prob = exp(-IMMOTILITY*dH);
				//if (prob>mt_random_double())
        		if (prob>(rand()/(double)RAND_MAX))
				{
            		pv[xt].ctag = stag; // a move is made
            		if(ttag) {csize[ttag-1]--;}
            		if(stag) {csize[stag-1]++;}
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
/* NO LONGER IN USE: non-local, for MPI want local alternative
// This function checks if retracting a voxel from a cell would cause it to split
// 	returns true if the retraction would split the cell (not allowed!)
//	returns false if the retraction is safe
BOOL splitcheckCCR(VOX* pv, int* csize, int xt, int ttag)
{
	BOOL split;
	int nbs[8],n,nb,prev,curr,in;
	int v, nrblue, nrgrey, startnb;
	int greys[csize[ttag-1]];
	short *CCAlabels = malloc(my_nvy * NVX *sizeof(short));
	int i, nrgrey0, g, nbsg[8];

	// get 8 neighbors of voxel xt
	nbs[0]=xt-1+NVX; nbs[1]=xt+NVX; nbs[2]=xt+1+NVX;
	nbs[7]=xt-1;                    nbs[3]=xt+1;
	nbs[6]=xt-1-NVX; nbs[5]=xt-NVX; nbs[4]=xt+1-NVX;

	// circle through neighbors and count how many transitions there are
	prev = pv[nbs[7]].ctag; in = 0;
	for(n=0;n<8;n++)
	{
		curr = pv[nbs[n]].ctag;
		if((prev!=ttag)&&(curr==ttag))
			in++;
		prev = curr;
	}

	split = FALSE;
	// if there is more than one transition then a split can occur
	if(in>1)
	{
		// CONNECTED COMPONENT ALGORITHM Rene-style (CCR)
    	// connected checking "label":
    	// 0: blue;    neighbors of retracted element
    	// 1: white;   undiscovered
    	// 2: grey;    discovered but not finished processing
    	// 3: black;   finished processing

		for(v=0;v<my_nvy * NVX;v++) 
		{
			CCAlabels[v] = 1; //set all voxel labels to undiscovered
		} 
		CCAlabels[xt] = 3; //voxel xt is processed

		// mark neighbors of retracted element
		nrblue = -1;
		for(n=0;n<8;n++)
		{
			nb = nbs[n];
			// if (nb < 0 || nb >= my_nvy * NVX){
			// 	printf("Rank %d; issue: sees value of nb = %d with nrblue = %d",myRank,nb,nrblue); fflush(stdout);
			// }
			// assert(nb >= 0 && nb < my_nvy * NVX);
			if(pv[nb].ctag==ttag)
			{
				CCAlabels[nb]=0; nrblue++;
				startnb = nb;
			}
		}

		CCAlabels[startnb]=2;
		nrgrey=1;
		greys[0]=startnb;

		while(nrgrey&&nrblue)
		{
			nrgrey0 = nrgrey;
			// make neighbors of discovered greys grey
			for(i=0;i<nrgrey0;i++)
			{
				g = greys[i];
				nbsg[0]=g-1+NVX; nbsg[1]=g+NVX; nbsg[2]=g+1+NVX;
				nbsg[7]=g-1;                    nbsg[3]=g+1;
				nbsg[6]=g-1-NVX; nbsg[5]=g-NVX; nbsg[4]=g+1-NVX;
				for(n=0;n<8;n++)
				{
					nb = nbsg[n];
					// if (nb < 0 || nb >= my_nvy * NVX){
					// 	printf("Rank %d; issue: sees value of nb = %d for g = %d",myRank,nb,g); fflush(stdout);
					// }
					// assert(nb >= 0 && nb < my_nvy * NVX);
					if((pv[nb].ctag==ttag)&&(CCAlabels[nb]<2))
					{
						if(CCAlabels[nb]==0) {nrblue--;}
						CCAlabels[nb]=2; nrgrey++; greys[nrgrey-1]=nb;
					}
				}
			}

			// make processed greys black
			for(i=0;i<nrgrey0;i++)
			{
				g = greys[i];
				CCAlabels[g]=3;
				greys[i]=greys[nrgrey-1]; nrgrey--;
			}

		}
		// if any blue voxels remain undiscovered, the cell would split
		if(nrblue) {split = TRUE;}

	}
	free(CCAlabels);
	return split;
}*/

// inspired by Durand and Guesnet (2016) https://doi.org/10.1016/j.cpc.2016.07.030
BOOL splitcheck_local(VOX* pv, int xt, int ttag)
{
    int x = xt % NVX;
    int y = xt / NVX;
	int nbs[8]; 

	// define neighbors
	nbs[0]=xt-1+NVX; nbs[1]=xt+NVX; nbs[2]=xt+1+NVX;
    nbs[7]=xt-1;                    nbs[3]=xt+1;
    nbs[6]=xt-1-NVX; nbs[5]=xt-NVX; nbs[4]=xt+1-NVX;

	// check for any Von Neumann neighbors
	int N = (pv[nbs[1]].ctag == ttag);
    int S = (pv[nbs[5]].ctag == ttag);
    int E = (pv[nbs[3]].ctag == ttag);
    int W = (pv[nbs[7]].ctag == ttag);

    int z = N + S + E + W;

    // trivial cases
    if (z <= 1) return FALSE;   // cannot split if only one pixel
    if (z == 4) return TRUE;    // removing center causes split for sure

    // check for diagonal neighbors
    int NE = (pv[nbs[2]].ctag == ttag);
    int NW = (pv[nbs[0]].ctag == ttag);
    int SE = (pv[nbs[4]].ctag == ttag);
    int SW = (pv[nbs[6]].ctag == ttag);

    if (z == 2)
    {
        // if opposite neighbors, would cause split
        if ((N && S) || (E && W))
            return TRUE;

        // if corner neighbors, would cuase split
        if (N && E && !NE) return TRUE;
        if (E && S && !SE) return TRUE;
        if (S && W && !SW) return TRUE;
        if (W && N && !NW) return TRUE;

        return FALSE;
    }

    if (z == 3)
    {
        // if missing one direction, check both diagonals bridging gap
        if (!N && (!SE || !SW)) return TRUE;
        if (!S && (!NE || !NW)) return TRUE;
        if (!E && (!NE || !SE)) return TRUE;
        if (!W && (!NW || !SW)) return TRUE;
        return FALSE;
    }

    return FALSE;
}


void exchange_cell_halos(VOX* pv, int NRc, int *NRC_loc, int **loc_cells)
{
	int c;
	int v, vx, vy;

	// Halo exchange
	// rank sends first two owned rows to topNbr, and botNbr receives into its bottom two rows (the halo)
	MPI_Sendrecv(&pv[2 * NVX], 2 * NVX, MPI_INT, topNbr, 0, &pv[(my_nvy-2) * NVX], 2 * NVX, MPI_INT, botNbr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// rank sends last two owned rows to botNbr, and topNbr receives into its top two rows (the halo)
	MPI_Sendrecv(&pv[(my_nvy-4) * NVX], 2 * NVX, MPI_INT, botNbr, 1, &pv[0], 2 * NVX, MPI_INT, topNbr, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// update loc_cells after halo exchange
	int *seen = calloc(NRc, sizeof(int)); // track whether or not cell ID has already been seen
	int loc_idx = 0;
    for (vy = 0; vy < my_nvy; vy++)	// loop through all local voxels
    for (vx = 0; vx < NVX; vx++)
    {
        v = vx + vy * NVX;
		c = pv[v].ctag;
        if (c > 0 && !seen[c-1]) {
            (*loc_cells)[loc_idx++] = pv[v].ctag; //store global cell ID in local list
			seen[c-1]=1;
        }
    }
	*NRC_loc = loc_idx; // save local cells after halo exchange
	free(seen);

	return;
}

void update_cell_volumes(VOX* pv, int* csize, int NRc)
{
	int v, vx, vy;

	memset(csize, 0, NRc*sizeof(int));
	for(vy=2; vy<my_nvy-2; vy++) //owned voxels
   	for(vx=0; vx<NVX; vx++)
   	{
   		v = vx + vy*NVX;
		if(pv[v].ctag > 0) { //if voxel v belongs to a cell (using global cell index), add it to that cell's volume
			csize[pv[v].ctag-1]++;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, csize, NRc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	// for debugging
	// for(int i =0; i<NRc; i++) //owned voxels
   	// {
	// 	printf("\nRank %d, sees csize value of %d",myRank,csize[i]);
	// }
	return;
}