# HPC Stencil Method Project on Orfeo

This project implements a hybrid MPI+OpenMP stencil method for solving the 2D Heat Equation, optimized for the Orfeo cluster.

## Project Structure
- `src/`: C source files.
- `include/`: Header files.
- `scripts/`: SLURM batch scripts for scaling tests.
- `Makefile`: Build configuration.

## How to Run on Orfeo

### 1. Connect to Orfeo
SSH into the Orfeo login node.

### 2. Build the Project
Load the necessary modules and run make:
```bash
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0
make
```

### 3. Run Scaling Tests
Submit the provided SLURM scripts to the scheduler:

#### OpenMP Scaling (Single Node)
```bash
sbatch scripts/run_omp_scaling.sh
```

#### MPI Strong Scaling (Multi-Node)
```bash
sbatch scripts/run_mpi_strong_scaling.sh
```

#### MPI Weak Scaling (Multi-Node)
```bash
sbatch scripts/run_mpi_weak_scaling.sh
```

## Results
The scripts will generate CSV files (`omp_results.csv`, `mpi_strong_results.csv`, `mpi_weak_results.csv`) containing the timing data. You can use these to plot your scaling curves.
