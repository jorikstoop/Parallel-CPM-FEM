// file: nodelabels.c
#include "functions.h"
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
void disp_to_nodes(NOD* pn, double* u)
{
	int n, nx, ny;
	int global_ny, owned;
	int cnt =0;

	for (ny = 0; ny < my_nny; ny++) //loop through nodes
	for (nx = 0; nx < NNX; nx++) {
		n = nx + ny * NNX;

		// check if node is owned by current rank
		global_ny = my_ymin + (ny - 2);
		owned = is_owned_node(global_ny);

        if (!owned) continue; //only consider owned nodes

		//assert(cnt < nrrdof_owned);

		pn[n].ux=.0; //reset to zero
		pn[n].uy=.0;
		if(!pn[n].restrictx)
			pn[n].ux=u[2*n];
		if(!pn[n].restricty)
			pn[n].uy=u[2*n+1];
	}
}

////////////////////////////////////////////////////////////////////////////////
void set_disp_of_prev_incr(NOD* pn, double* u)
{
	int n, nx, ny;
	int global_ny, owned;
	int cnt =0;

	for (ny = 0; ny < my_nny; ny++) // loop through nodes
	for (nx = 0; nx < NNX; nx++) {
		n = nx + ny * NNX;

		// check if node is owned by current rank
		global_ny = my_ymin + (ny - 2);
		owned = is_owned_node(global_ny);

        if (!owned) continue; //only consider owned nodes
		
		if(!pn[n].restrictx)
			u[2*n]=pn[n].ux;
		if(!pn[n].restricty)
			u[2*n+1]=pn[n].uy;
	}
}

////////////////////////////////////////////////////////////////////////////////
void place_node_forces_in_f(NOD* pn, double* f)
{
	int n, nx, ny;
	int global_ny, owned;
	int cnt =0;

	for (ny = 0; ny < my_nny; ny++) //loop through nodes
	for (nx = 0; nx < NNX; nx++) {

		n = nx + ny * NNX; //loc node idx

		// check if node is owned by current rank
		global_ny = my_ymin + (ny - 2);
		owned = is_owned_node(global_ny);

        if (!owned) continue; //only consider owned nodes
        		
		if(!pn[n].restrictx)
			f[2*n]=pn[n].fx;
		if(!pn[n].restricty)
			f[2*n+1]=pn[n].fy;
	}

}

