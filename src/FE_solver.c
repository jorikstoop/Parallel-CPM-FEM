#include "functions.h"
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
//new Kdotx with csr
void calc_Kdotx_csr(CSRMatrix* A, double* x, double* b)
{
    int r, j;

    #pragma omp parallel for private(r,j)
    for (r = 0; r < A->n; r++) {

        int n = r / 2;
        int ny = n / NNX;
        int global_ny = my_ymin + (ny - 2);

        if (!is_owned_node(global_ny)) continue;

        double sum = 0.0;

        for (j = A->rowptr[r]; j < A->rowptr[r+1]; j++) {
            sum += A->vals[j] * x[A->colind[j]];
        }

        b[r] = sum;
    }
}
//new PCGsolve with csr format
////////////////////////////////////////////////////////////////////////////////
void solvePCG_csr(CSRMatrix* A, double* u, double* f, NOD* pn)
{
    int i, k, iter;
    double *ui, *ri, *invC, *zi, *pi, *qi;
    double rhoi, rhoinew, initrho;
    double beti, alfi, pq;

    ui   = calloc(my_ndof, sizeof(double));
    ri   = calloc(my_ndof, sizeof(double));
    invC = calloc(my_ndof, sizeof(double));
    zi   = calloc(my_ndof, sizeof(double));
    pi   = calloc(my_ndof, sizeof(double));
    qi   = calloc(my_ndof, sizeof(double));

    // create array to store the indices of owned, un-restricted dofs to speed up subsequent for loops
    int n, nx, ny, global_ny, owned, base;
    int *owned_dofs = malloc(my_ndof * sizeof(int));
    int n_owned = 0;
    for (ny = 0; ny < my_nny; ny++)
    for (nx = 0; nx < NNX; nx++)
    {
        n = nx + ny * NNX;

        // check if node is owned
		global_ny = my_ymin + (ny - 2);
		owned = is_owned_node(global_ny);
        if (!owned) continue;

        base = 2*n;
        if (!pn[n].restrictx) owned_dofs[n_owned++] = base;
        if (!pn[n].restricty) owned_dofs[n_owned++] = base + 1;
    }

    //Diagonal preconditioner (extract from CSR)
    #pragma omp parallel for
    for (i = 0; i < my_ndof; i++) { //initialize over all local dofs including halos
        ui[i] = u[i]; 
        invC[i] = (A->diag[i] != 0.0) ? 1.0 / A->diag[i] : 0.0;
    }

    exchange_vector_halos(ui); //is this needed since we memset u to zero every iter?
    calc_Kdotx_csr(A, ui, qi);

    rhoinew = 0.0;
    initrho = 0.0;

    #pragma omp parallel for private(k,i) reduction(+:rhoinew, initrho)
    for (k = 0; k < n_owned; k++) {
        i = owned_dofs[k];

        ri[i] = f[i] - qi[i];
        zi[i] = invC[i] * ri[i];
        pi[i] = zi[i];
        rhoinew += ri[i] * zi[i];
        initrho += invC[i] * f[i] * f[i];
    }
    // reduce global values
    MPI_Allreduce(MPI_IN_PLACE, &rhoinew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &initrho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for (iter = 0; rhoinew > ACCURACY * initrho; iter++) {

        rhoi = rhoinew;

        exchange_vector_halos(pi); // exchange halos first
        calc_Kdotx_csr(A, pi, qi);

        pq = 0.0;
        #pragma omp parallel for private(k,i) reduction(+:pq)
        for (k = 0; k < n_owned; k++) {
            i = owned_dofs[k];

            pq += pi[i] * qi[i];
        }
        //pq is global
        MPI_Allreduce(MPI_IN_PLACE, &pq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        alfi = rhoi / pq;

        #pragma omp parallel for private(k,i)
        for (k = 0; k < n_owned; k++) {
            i = owned_dofs[k];

            ui[i] += alfi * pi[i];
            ri[i] -= alfi * qi[i];
        }

        rhoinew = 0.0;
        #pragma omp parallel for private(k,i) reduction(+:rhoinew)
        for (k = 0; k < n_owned; k++) {
            i = owned_dofs[k];

            zi[i] = invC[i] * ri[i];
            rhoinew += ri[i] * zi[i];
        }

        //rhoinew is global
        MPI_Allreduce(MPI_IN_PLACE, &rhoinew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        beti = rhoinew / rhoi;

        #pragma omp parallel for private(k,i)
        for (k = 0; k < n_owned; k++) {
            i = owned_dofs[k];

            pi[i] = zi[i] + beti * pi[i];
        }
        if (myRank == 0 && iter % 10 == 0)
            printf("\ni %4d, rhoinew/initrho=%18.11lf", iter, rhoinew / initrho);
    }

    if (myRank == 0) printf("\n Stop iterating at iter %d", iter);

    #pragma omp parallel for private(k,i)
        for (k = 0; k < n_owned; k++) {
            i = owned_dofs[k];

        u[i] = ui[i];
    }

    free(ui);
    free(ri);
    free(invC);
    free(zi);
    free(pi);
    free(qi);
    free(owned_dofs);
}

void exchange_vector_halos(double *v)
{
    double temp = MPI_Wtime();

    // note: only send and recv one row of nodes bc only 1 layer halo needed for FEM
    int count = 2 * NNX; // 1 row of nodes = NNX nodes × 2 dof


    int last_owned_ny, inner_bot_halo_ny;
    if (myRank == size - 1){
        last_owned_ny = my_nny - 3;
        inner_bot_halo_ny = my_nny - 2;
    }else{
        last_owned_ny = my_nny - 4;
        inner_bot_halo_ny = my_nny - 3;
    }
	//printf("Rank %d, last_owned_ny =  %d \n", myRank, last_owned_ny);

    // send top owned row of nodes to topNbr, who receives into their bottom halo
    // (note that because we assign ownership of shared nodes to upper neighbor, the first row of the bottom halo is actually my_nny-3)
    MPI_Sendrecv(&v[2*(2*NNX)], count, MPI_DOUBLE, topNbr, 0,
                 &v[2*(inner_bot_halo_ny*NNX)], count, MPI_DOUBLE, botNbr, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // send bottom owned row to botNbr, who receives into their top halo 
    // (note that because we assign ownership of shared nodes to upper neighbor, the last owned row is actually my_nny-4)
    MPI_Sendrecv(&v[2*(last_owned_ny*NNX)], count, MPI_DOUBLE, botNbr, 1,
                 &v[2*(1*NNX)], count, MPI_DOUBLE, topNbr, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    time_fem_comm += MPI_Wtime() - temp;
}

void exchange_node_halos(NOD* pn)
{
    // note: only send and recv one row of nodes bc only 1 layer halo needed for FEM
    int count = NNX;

    int last_owned_ny, inner_bot_halo_ny;
    if (myRank == size - 1){
        last_owned_ny = my_nny - 3;
        inner_bot_halo_ny = my_nny - 2;
    }else{
        last_owned_ny = my_nny - 4;
        inner_bot_halo_ny = my_nny - 3;
    }

    // send top owned row of nodes to topNbr, and receive into bottom halo 
    // (note that because we assign ownership of shared nodes to upper neighbor, the first row of the bottom halo is actually my_nny-3)
    MPI_Sendrecv(&pn[2*NNX], count*sizeof(NOD), MPI_BYTE, topNbr, 0,
                 &pn[inner_bot_halo_ny*NNX], count*sizeof(NOD), MPI_BYTE, botNbr, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // send bottom owned row to botNbr and receive into top halo
    // (note that because we assign ownership of shared nodes to upper neighbor, the last owned row is actually my_nny-4)
    MPI_Sendrecv(&pn[last_owned_ny*NNX], count*sizeof(NOD), MPI_BYTE, botNbr, 1,
                 &pn[1*NNX], count*sizeof(NOD), MPI_BYTE, topNbr, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

/*
void calc_Kdotx(int* kcol, double* kval, double* diag, double* x, double* b, int nrrdof)
{
    int r, a, lim, tid;
    int nthreads;
    double **b_private;

    nthreads = omp_get_max_threads();

    b_private = (double **)malloc(nthreads * sizeof(double *));
    if (b_private == NULL) {
        printf("Error allocating b_private\n");
        exit(1);
    }

    for (tid = 0; tid < nthreads; tid++) {
        b_private[tid] = (double *)calloc(nrrdof, sizeof(double));
        if (b_private[tid] == NULL) {
            printf("Error allocating b_private[%d]\n", tid);
            exit(1);
        }
    }

    #pragma omp parallel private(tid, r, a, lim)
    {
        tid = omp_get_thread_num();

        #pragma omp for
        for (r = 0; r < nrrdof; r++) {
            b_private[tid][r] = diag[r] * x[r];
        }

        #pragma omp for
        for (r = 0; r < nrrdof; r++) {
            lim = 10 * r + kcol[10 * r];
            for (a = 10 * r + 1; a < lim; a++) {
                int c = kcol[a];
                double val = kval[a];
                b_private[tid][r] += val * x[c];
                b_private[tid][c] += val * x[r];
            }
        }
    }

    #pragma omp parallel for private(tid)
    for (r = 0; r < nrrdof; r++) {
        b[r] = 0.0;
        for (tid = 0; tid < nthreads; tid++) {
            b[r] += b_private[tid][r];
        }
    }

    for (tid = 0; tid < nthreads; tid++) {
        free(b_private[tid]);
    }
    free(b_private);
}
*/
/*
void solvePCG(int* kcol, double* kval, double* u, double* f, int nrrdof)
{
    int i, iter;
    double *ui, *ri, *diag, *invC, *zi, *pi, *qi;
    double rhoi, rhoinew, initrho;
    double beti, alfi, pq;

    ui   = malloc(nrrdof * sizeof(double));
    ri   = malloc(nrrdof * sizeof(double));
    diag = malloc(nrrdof * sizeof(double));
    invC = malloc(nrrdof * sizeof(double));
    zi   = malloc(nrrdof * sizeof(double));
    pi   = malloc(nrrdof * sizeof(double));
    qi   = malloc(nrrdof * sizeof(double));

    #pragma omp parallel for
    for (i = 0; i < nrrdof; i++) {
        ui[i] = u[i];
        diag[i] = kval[10 * i];
        if (diag[i] != 0.0)
            invC[i] = 1.0 / diag[i];
        else
            invC[i] = 0.0;
    }

    calc_Kdotx(kcol, kval, diag, ui, qi, nrrdof);

    rhoinew = 0.0;
    initrho = 0.0;
    #pragma omp parallel for reduction(+:rhoinew, initrho)
    for (i = 0; i < nrrdof; i++) {
        ri[i] = f[i] - qi[i];
        zi[i] = invC[i] * ri[i];
        pi[i] = zi[i];
        rhoinew += ri[i] * zi[i];
        initrho += invC[i] * f[i] * f[i];
    }

    for (iter = 0; rhoinew > ACCURACY * initrho; iter++) {
        rhoi = rhoinew;

        calc_Kdotx(kcol, kval, diag, pi, qi, nrrdof);

        pq = 0.0;
        #pragma omp parallel for reduction(+:pq)
        for (i = 0; i < nrrdof; i++) {
            pq += pi[i] * qi[i];
        }

        alfi = rhoi / pq;

        #pragma omp parallel for
        for (i = 0; i < nrrdof; i++) {
            ui[i] += alfi * pi[i];
            ri[i] -= alfi * qi[i];
        }

        rhoinew = 0.0;
        #pragma omp parallel for reduction(+:rhoinew)
        for (i = 0; i < nrrdof; i++) {
            zi[i] = invC[i] * ri[i];
            rhoinew += ri[i] * zi[i];
        }

        beti = rhoinew / rhoi;

        #pragma omp parallel for
        for (i = 0; i < nrrdof; i++) {
            pi[i] = zi[i] + beti * pi[i];
        }

        if (iter % 10 == 0)
            printf("\ni %4d, rhoinew/initrho=%18.11lf", iter, rhoinew / initrho);
    }

    printf("\n Stop iterating at iter %d", iter);

    #pragma omp parallel for
    for (i = 0; i < nrrdof; i++) {
        u[i] = ui[i];
    }

    free(ui);
    free(ri);
    free(diag);
    free(invC);
    free(zi);
    free(pi);
    free(qi);
}
	*/