#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=EPYC
#SBATCH --account=dssc
#SBATCH --time=00:05:00
#SBATCH --job-name=test_stencil
#SBATCH --output=test_%j.out
#SBATCH --error=test_%j.err

module purge
module load openMPI/5.0.5

echo "========================================="
echo "ORFEO Test Run"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "========================================="

# Check if executable exists
if [ ! -f ./stencil_hybrid ]; then
    echo "ERROR: Executable ./stencil_hybrid not found!"
    echo ""
    echo "Compile with:"
    echo "  module load openMPI/5.0.5"
    echo "  mpicc -fopenmp -O3 -march=native -o stencil_hybrid stencil_parallel.c -lm"
    exit 1
fi

# OpenMP settings
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_DISPLAY_ENV=TRUE

echo ""
echo "Test Configuration:"
echo "  - 1 MPI task"
echo "  - ${OMP_NUM_THREADS} OpenMP threads"
echo "  - Grid: 500x500"
echo "  - Iterations: 100"
echo ""

srun --ntasks=1 --cpus-per-task=8 ./stencil_hybrid -x 500 -y 500 -n 100

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "Test completed successfully!"
    echo "========================================="
else
    echo ""
    echo "========================================="
    echo "Test FAILED!"
    echo "========================================="
    exit 1
fi
