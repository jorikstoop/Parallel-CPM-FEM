#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
/*this handles matrix stoage and formatting :
conversts 
`kcol (indices)
kval (values)`

into 

`rowptr (size nrrdof+1)
colind (size nnz)
vals   (size nnz)`
*/
CSRMatrix convert_to_csr(int* kcol, double* kval)
{
    int r, c, a, lim;
    int *row_counts;
    CSRMatrix A;

    //count nnz per row
    row_counts = calloc(my_ndof, sizeof(int));

    for (r = 0; r < my_ndof; r++) {
        row_counts[r]++; //diagonal

        lim = 10*r + kcol[10*r];
        for (a = 10*r + 1; a < lim; a++) {
            c = kcol[a]; // could be owned or halo
            row_counts[r]++;
            row_counts[c]++;
        }
    }

    //allocate CSR arrays
    A.n = my_ndof;
    A.rowptr = malloc((my_ndof + 1) * sizeof(int));

    A.rowptr[0] = 0;
    for (r = 0; r < my_ndof; r++) {
        A.rowptr[r+1] = A.rowptr[r] + row_counts[r];
    }

    A.nnz = A.rowptr[my_ndof];

    A.colind = malloc(A.nnz * sizeof(int));
    A.vals   = malloc(A.nnz * sizeof(double));

    //allocate diagonal storage
    A.diag = malloc(my_ndof * sizeof(double));

    //3: reset counters
    for (r = 0; r < my_ndof; r++) {
        row_counts[r] = 0;
    }

    //4: fill CSR arrays
    for (r = 0; r < my_ndof; r++) {

        //diagonal
        int idx = A.rowptr[r] + row_counts[r];
        A.colind[idx] = r;
        A.vals[idx]   = kval[10*r];
        A.diag[r]     = kval[10*r];   //store diagobal
        row_counts[r]++;

        lim = 10*r + kcol[10*r];

        for (a = 10*r + 1; a < lim; a++) {
            c = kcol[a];
            double val = kval[a];

            //(r, c)
            idx = A.rowptr[r] + row_counts[r];
            A.colind[idx] = c;
            A.vals[idx]   = val;
            row_counts[r]++;

            //(c, r)
            idx = A.rowptr[c] + row_counts[c];
            A.colind[idx] = r;
            A.vals[idx]   = val;
            row_counts[c]++; //this DOUBLES STORAGE
        }
    }

    free(row_counts);

    return A;
}