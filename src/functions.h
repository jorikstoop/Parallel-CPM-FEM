// file: functions.h
#include "def.h"
#include "structures.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>

// global variables
extern int myRank, size;
extern int my_ymin, my_ymax;
extern int my_nvy, my_nny, my_ndof;
extern int my_nvy_owned, my_nny_owned, my_ndof_owned;
extern int topNbr, botNbr;
extern double time_cell_comm;
extern double time_fem_comm;

// init.c
VOX*		init_voxels(void);
NOD*		init_nodes(void);
uint32_t    hash_vox(int gvx, int gvy, int seed);
int         init_cells(VOX* pv, int *NRc_loc, int **loc_cells);
void 		set_forces(NOD* pn);
void		set_restrictions(NOD* pn);
int         is_owned_node(int global_ny);

// cellmoves.c
void 		CPM_moves(VOX* pv, NOD* pn, int* csize);
BOOL 		splitcheckCCR(VOX* pv,  int* csize, int xt, int ttag);
BOOL        splitcheck_local(VOX* pv, int xt, int ttag);
void        exchange_cell_halos(VOX* pv, int NRc, int *NRc_loc, int **loc_cells);
void        update_cell_volumes(VOX* pv, int* csize, int NRc);

// CPM_dH.c
double 		calcdH(VOX* pv, NOD* pn, int* csize, int xt, int xs, int pick, int ttag, int stag);
double 		calcdHcontact(VOX* pv, int xt, int ttag, int stag);
double 		contactenergy(int tag1, int tag2);
double 		calcdHvol(int* csize, int ttag, int stag);
double 		calcdHstrain(NOD* pn, int xt, int xs, int pick, int ttag, int stag);
double 		sige(double L);


// cellforces.c
void 		cell_forces(VOX* pv, NOD* pn, int nthreads, int NRc, int NRc_loc, int* loc_cells);


// FE_local.c
double** 	set_klocal(void);
void 		material_matrix(double *pD);
void 		set_matrix_B(double *pB, double x, double y);
void 		get_estrains(NOD* pn, int vx, int vy, double* estrains);
void 		get_estress(int e, double* estrains, double* estress);
void 		get_princs(double* str, double* L1, double* L2, double* v1, double* v2, BOOL strain);

// FE_assembly.c
//int** 		set_topology(void);
void		assembly(int* kcol, double* kval, double** klocal, VOX* pv);
int 		arrange_dofpos(int* dofpos, NOD* pn, int* nrrdof_owned);
void 		reduce_K(int* kcol, double* kval, int* dofpos);

// FE_nodes2dofs.c
void 		disp_to_nodes(NOD* pn, double* u);
void 		set_disp_of_prev_incr(NOD* pn, double* u);
void 		place_node_forces_in_f(NOD* pn, double* f);

// FE_csr.c
// please see FE_csr.c
typedef struct {
    int n;
    int nnz;
    int *rowptr;
    int *colind;
    double *vals;
    double *diag; 
} CSRMatrix;

CSRMatrix convert_to_csr(int* kcol, double* kval);


// FE_solver.c
void        calc_Kdotx_csr(CSRMatrix* A, double* x, double* b);
void 		calc_Kdotx(int* kcol, double* kval, double* diag, double* x, double* b, int nrrdof);
void 		solvePCG(int* kcol, double* kval, double* u, double* f, int nrrdof);
void        solvePCG_csr(CSRMatrix* A, double* u, double* f, NOD* pn);
void        exchange_vector_halos(double *v);
void        exchange_node_halos(NOD* pn);

// write.c
void   		write_increment(int increment);
void 		write_cells(VOX* pv, int increment);
void 		write_pstrain(VOX* pv, NOD* pn, int increment);
void 		write_pstress(VOX* pv, NOD* pn, int increment);
void 		write_forces(NOD* pn, int increment);
void 		write_disps(NOD* pn, int increment);
void        write_timing(double time_setup, double time_simloop,double time_output, double time_tracforce,
                        double time_fem, double time_pcg, double time_fem_comm,
                        double time_cpm, double time_cell_comm, int increment);

// read.c
int   		read_increment(void);
int 		read_cells(VOX* pv, int increment);

// mylib.c
void 		myitostr(int n, char s[]);
void 		myreverse(char s[]);
unsigned 	mystrlen(const char *s);

// mt.c
void 		mt_init();
unsigned long mt_random();





