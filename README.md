# Parallel-CPM-FEM

Hybrid OpenMP-MPI parallel cellular Potts and finite element code for simulating endothelial cell behavior. This project was completed as part of the BME 520L course at Duke University.

Authors: Jorik Stoop, Guinevere Ferreira

### Attribution
This repository extends the simulation code from:

van Oers RFM, Rens EG, LaValley DJ, Reinhart-King CA, Merks RMH (2014).
“Mechanical Cell-Matrix Feedback Explains Pairwise and Collective Endothelial
Cell Behavior In Vitro.” *PLOS Computational Biology* 10(8): e1003774.
https://doi.org/10.1371/journal.pcbi.1003774

The original serial C and MATLAB code can be found [here](https://doi.org/10.1371/journal.pcbi.1003774.s006), and is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)

Changes made in this repository include parallel implementation using MPI and OpenMP, as well as Python postprocessing tools for visualization and analysis of simulation output. These changes were made by Jorik Stoop and Guinevere Ferreira.


## Compiling and Running the Code

### 1. Clone repo
- First clone the repository from GitHub onto your system
- Then use `cd Parallel-CPM-FEM` to enter the project directory

### 2. Setup environment
- Make sure you have GNU Make and a GNU C Compiler installed on your system
- Most HPC systems will have these automatically installed. Use `which cc` and `cc --version` to check
- Make sure you have MPI installed. Check with `which mpicc`
  - On the Duke Compute Cluster use `module load OpenMPI/4.1.6` to load MPI

### 3. Compile code
- Navigate to the `./src` folder with `cd src`
- Type `make` into the terminal to compile
- NOTE: there is a batch script which can be run from the main project directory using `source recompile.sh`

### 4. Run the parallel code
- Type `mpirun -np r ./src/cpmfem` (where r is the number of ranks) from the main project directory (not the src directory)

Example output:
```
Number of omp threads on rank 0 set to 2
Number of omp threads on rank 1 set to 2
Rank 0: global rows [0, 150), my_nvy=154, my_nny=155
Rank 1: global rows [150, 300), my_nvy=154, my_nny=155
Rank 0: sees 476 total cells in the simulation, with 241 local cells (including halo) 
Rank 1: sees 476 total cells in the simulation, with 240 local cells (including halo) 

ASSEMBLY COMPLETED
Total nrrdof owned = 178802
START INCREMENT 0
i    0, rhoinew/initrho=     0.06988641335
i   10, rhoinew/initrho=     0.00020242563
i   20, rhoinew/initrho=     0.00003214867
i   30, rhoinew/initrho=     0.00001180028
 Stop iterating at iter 33
```
The output log will show which rank owns which rows of the simulated domain.

For larger jobs (e.g. large NVX value to increase simulation domain size) this code can be run on HPC systems with multiple nodes. Here is an example SLURM run script:
```
#!/bin/bash
#SBATCH --nodes=8               # Total # of nodes
#SBATCH --ntasks-per-node=16    # Total # of MPI tasks per node
#SBATCH --cpus-per-task=2       # cpu-cores per task
#SBATCH --time=1:00:00          # Total run time limit (hh:mm:ss)
#SBATCH -J cpm-fem              # Job name
#SBATCH -o cpm-fem.o%j          # Name of stdout output file
#SBATCH -e cpm-fem.e%j          # Name of stderr error file
#SBATCH -p <partition_name>     # change this based on system

module load OpenMPI/4.1.6   # change this based on system
source recompile.sh

# Set OpenMP thread count
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Launch hybrid MPI + OpenMP code
mpirun -np 128 ./src/cpmfem
```


## Visualization of endothelial network formation

### 1. Setup Python environment
- we recommend using a virtual environment
- load the required modules from the requirements.txt file with `pip3 install -r requirements.txt`

### 2. Make animation
- enter the `./postprocess` directory with `cd postprocess`
- use `python3 make_movie.py --stop <incr> --nex <size>` to start the code
  - stop increment designates the end time step for the animation
  - NEX is the domain side length 
  - This script takes in the strain and cell position values from the `.out` files found in the `output` folder to create a mp4 movie file 

NOTE: there are additional command-line arguments for the `make_movie.py` code listed below:

- "--path" or "-p", default="../output", Directory with ctags#.out and pstrain#.out files
- "--start", default=0, start increment (inclusive)
- "--stop", default=3000, Stop increment (exclusive)
- "--step", default=5, Stride between increments
- "--out" or "-o", default="simmovie.mp4", Output movie filename
- "--fps", default=25, Frames per second
- "--nex", default=None, NEX
- "--ney", default=None, NEY (optional)

## Run Stiffness Study

### Comparing parallel implementation results to published data

In the source publication, the authors perform a series of tests observing single cell behavior on ECM with varying stiffness. The authors measured the cell area and length after 100 time steps on substrates of stiffness 0.5, 1, 2, 4, 8, 10, 12, 14, 16, and 32 kPa.

To verify our parallel implementation, we repeat the study.

### Simulation setup
To setup these simulations, first modify the following parameters in `./src/def.h`:
- set the domain size to 100 with `#define NVX 100`
- set the number of steps to 100 with `#define NRINC 100`
- set the stiffness with `#define YOUNGS 10E3 // [Pa]`. Note that units here are in Pa, so 10E3 is 10kPa

Then change the `./src/init.c` function to only add one cell at the center of the domain:
- comment out the line `if (r01 < .25 / TARGETVOLUME) // place cells with 25% of max density` in the init_cells() function.
- uncomment the line `if((vx==NVX/2)&&(vy==NVY/2)) // to place cell at exact global center` in the init_cells() function.

Once those changes are made, you can run the stiffness parameter sweep with `source run_stiffness_study.sh`. This command runs all ten simulations in a row. To run in parallel, just update "MPI_RANKS=r" and "export OMP_NUM_THREADS=n" to the desired number of ranks and threads per rank before starting the parameter sweep.

Note that `N_REPS` is a variable that can be changed to alter the number of replicate runs performed at each stiffness value. The source paper uses N_REPS=100. Here we change the random SEED to be equal to N_REPS to introduce stochasticity across runs.

### Analyze result
After the simulations are complete, cell area & length can be measured using the `measure_cell_values.py` file in the postprocessing directory. 

This python script parses the output directories, calculates cell area (num pixels) and cell length (max dist between pixels) and then outputs two graphs: one for cell area vs stiffness and one for cell length vs stiffness. The resulting graphs should show an increase in cell area with increased stiffness, and a peak cell length on substrate stiffness of 12-14 kPa.


## Profiling
The code prints a timing breakdown to screen when the simulation is finished.

For additional analysis, the `profiling.py` script can be used by running `python3 profiling.py` from the postprocessing directory. This script parses the .csv timing files saved in the output directory and makes two graphs: a stacked bar graph showing the overall runtime, and a horizontal stacked bar showing the percentage breakdown of the run time.