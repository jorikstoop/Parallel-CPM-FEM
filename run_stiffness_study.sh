#!/bin/bash

# set num threads and MPI ranks
MPI_RANKS=4
export OMP_NUM_THREADS=1

# number of replicates
N_REP=3

# stiffness values (kPa)
STIFFNESSES=(0.5 1 2 4 8 10 12 14 16 32)

for KPA in "${STIFFNESSES[@]}"
do
    echo "=============================="
    echo "Running for ${KPA} kPa"
    echo "=============================="

    # convert kPa to Pa
    YOUNGS=$(awk "BEGIN {print ${KPA} * 1000}")

    echo "Setting YOUNGS = ${YOUNGS} Pa"

    # ---- modify def.h ----
    sed -i "s/#define YOUNGS .*/#define YOUNGS ${YOUNGS} \/\/ [Pa]/" src/def.h

    # run replicates
    for ((rep=1; rep<=N_REP; rep++))
    do
        echo "  -> Replicate ${rep}"

        # set SEED = replicate number
        sed -i "s/#define SEED.*/#define SEED ${rep}/" src/def.h

        # recompile
        source recompile.sh

        # clean old outputs
        rm -f output/*.out

        # run simulation
        mpirun -np $MPI_RANKS ./src/cpmfem

        # make directory
        OUTDIR="output/run_${KPA}_kPa_rep${rep}"
        mkdir -p "$OUTDIR"

        # move outputs
        mv output/*.out "$OUTDIR"/
    done
    echo ""
done

echo "All simulations completed."