#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=EPYC
#SBATCH --account=dssc
#SBATCH --time=00:30:00
#SBATCH --job-name=omp_scaling
#SBATCH --output=omp_scaling_%j.out
#SBATCH --error=omp_scaling_%j.err

# Load modules
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

EXEC="./stencil_hybrid"
ARGS="-x 2000 -y 2000 -n 1000"

# Thread counts for scaling study (1 to 128 for EPYC)
THREAD_COUNTS=(1 2 4 8 16 32 64 128)

# Ensure the results directory exists
mkdir -p results

echo "Threads,TotalTime,CommTime,CompTime" > results/omp_results.csv

for t in "${THREAD_COUNTS[@]}"; do
    export OMP_NUM_THREADS=$t
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo "========================================="
    echo "Running with $t OpenMP threads..."
    echo "========================================="
    
    # Run with single MPI task
    srun --ntasks=1 --cpus-per-task=$t $EXEC $ARGS > out_omp_$t.txt 2>&1
    
    # Extract timing information
    total=$(grep "Total time:" out_omp_$t.txt | awk '{print $3}')
    comm=$(grep "Communication time:" out_omp_$t.txt | awk '{print $3}')
    comp=$(grep "Computation time:" out_omp_$t.txt | awk '{print $3}')
    
    echo "$t,$total,$comm,$comp" >> results/omp_results.csv
    echo "Completed: $t threads - Total: ${total}s, Comm: ${comm}s, Comp: ${comp}s"
done

echo ""
echo "========================================="
echo "OpenMP scaling study complete!"
echo "Results saved to: results/omp_results.csv"
echo "========================================="