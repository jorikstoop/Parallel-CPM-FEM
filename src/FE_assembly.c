#include "functions.h"
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
/*
int** set_topology(void)
{
	int **top;
	int v, vx, vy;
	int n00, n10, n11, n01;

	top = calloc(NV,sizeof(int*));
	for(v=0;v<NV;v++)
		top[v] = calloc(8,sizeof(int));

	// set topology
   	for(vy=0; vy<NVY; vy++)
   	for(vx=0; vx<NVX; vx++)
   	{
   		v = vx + vy*NVX;

		// determine corner node numbers of this voxel
		n00 = (vx  ) + (vy  )*NNX;
		n10 = (vx+1) + (vy  )*NNX;
		n11 = (vx+1) + (vy+1)*NNX;
		n01 = (vx  ) + (vy+1)*NNX;

		top[v][0] = 2*n00;
		top[v][1] = 2*n00+1;
   		top[v][2] = 2*n10;
		top[v][3] = 2*n10+1;
		top[v][4] = 2*n11;
		top[v][5] = 2*n11+1;
		top[v][6] = 2*n01;
		top[v][7] = 2*n01+1;
	}

	return top;
}
*/

////////////////////////////////////////////////////////////////////////////////
void assembly(int* kcol, double* kval, double** klocal, VOX* pv)
{
	int v, vx, vy;
	int vy_global;
	double value;
	int il, jl, ig, jg; // row i & column j in klocal and K
	int d, a, b, lim;
	BOOL alreadynonzero;
	int n00, n10, n11, n01;
	int topv[8];
	double Ef; // multiply klocal with Ef depending on local E

	for(d=0;d<my_ndof;d++)
	{
		kcol[10*d] = 1;
		kval[10*d] = .0;
	}
 	
	// loop through voxels that include the surrounding 1-layer of node halos for FEM
	int vox_end = (myRank == size - 1) ? my_nvy-1 : my_nvy-2;
	for(vy=1; vy<vox_end; vy++)
	for(vx=0; vx<NVX; vx++)
	{
		Ef = 1;//(vx+1)/(double)NVX+.5;

		// determine corner node numbers of this element
		n00 = (vx  ) + (vy  )*NNX;
		n10 = (vx+1) + (vy  )*NNX;
		n11 = (vx+1) + (vy+1)*NNX;
		n01 = (vx  ) + (vy+1)*NNX;

		topv[0] = 2*n00;
		topv[1] = 2*n00+1;
   		topv[2] = 2*n10;
		topv[3] = 2*n10+1;
		topv[4] = 2*n11;
		topv[5] = 2*n11+1;
		topv[6] = 2*n01;
		topv[7] = 2*n01+1;


		// place klocal in K matrix
		for(il=0;il<8;il++) // go through rows in klocal
		for(jl=0;jl<8;jl++) // go through columns in klocal
		{
			value = Ef*klocal[il][jl];
			ig = topv[il]; // row in K
			jg = topv[jl]; // column in K

			if(jg==ig) // if on diagonal
				kval[10*ig] += value;

			if(jg>ig) // if right of diagonal
			{
				lim = 10*ig+kcol[10*ig];

				// check if there was already a nonzero on K(ig,jg)
				alreadynonzero = FALSE;
				for(a=10*ig+1;a<lim;a++) // go over ig-th row of K
				{
					if(kcol[a]==jg) // if storage for jg-th column of K
					{
						alreadynonzero = TRUE;
						b = a;
					}
				}

				if(alreadynonzero) // if already a nonzero on K(ig,jg)
					kval[b] += value; // add klocal(il,jl)

				else // if nothing on K(ig,jg)
				{
					b = lim;
					kcol[b] = jg; // make storage for jg-th column of K
					kval[b] = value; // add klocal(il,jl)
					kcol[10*ig]++;
				}
			} //endfor go through klocal
		} //endif relevant element
	} // endfor go though elements
	if (myRank == 0) printf("\nASSEMBLY COMPLETED");
}

////////////////////////////////////////////////////////////////////////////////
int arrange_dofpos(int* dofpos, NOD* pn, int* nrrdof_owned)
{
    int n, nx, ny;
	int global_ny, owned;
    int cnt_owned=0, cnt_tot=0;

	//fill dofpos based on spatial coords (instead of cnt)
    for(ny=0; ny<my_nny; ny++) // loop through all local nodes
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

		//check if node is owned by rank
		global_ny = my_ymin + (ny - 2);
		owned = is_owned_node(global_ny);

        // x DOF
        if (pn[n].restrictx){
            dofpos[2*n] = -1; //set restricted nodes to -1
		}
		else{
            dofpos[2*n] = 2*n; //set owned or halo node based on position n
         	if (owned) cnt_owned++;
			cnt_tot++;
		}

        // y DOF
        if (pn[n].restricty){
            dofpos[2*n+1] = -1;
		}	
        else{
			dofpos[2*n+1] = 2*n+1; //owned or halo node
         	if (owned) cnt_owned++;
			cnt_tot++;
		}
    }
	*nrrdof_owned = cnt_owned;

    return cnt_tot;
}


////////////////////////////////////////////////////////////////////////////////
void reduce_K(int* kcol, double* kval, int* dofpos)
{
	int ro, co; // row and column in K
	int a, lim, shift;

	for(ro=0;ro<my_ndof;ro++) //loop over all local DOFs
	{
		if(dofpos[ro]>-1) // exclude restricted rows
		{
			lim = 10*ro+kcol[10*ro];

			// change column indices for this row:
			for(a=10*ro+1;a<lim;a++)
			{
				co = kcol[a]; // old column index
				kcol[a] = dofpos[co];; // give new column index (can be owned, halo, or restricted)
			}
			// remove columns with -1 index (restricted)
			shift = 0;
			for(a=10*ro+1;a<lim;a++)
			{
				if(kcol[a]==-1){
					shift++;
				} else{
					kcol[a-shift] = kcol[a];
					kval[a-shift] = kval[a];
				}
			}
			kcol[10*ro] = kcol[10*ro]-shift;
		}
	}
}