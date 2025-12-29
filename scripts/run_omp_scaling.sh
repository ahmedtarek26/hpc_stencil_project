#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=EPYC
#SBATCH --time=00:30:00
#SBATCH --job-name=omp_scaling
#SBATCH --output=omp_scaling_%j.out
#SBATCH --error=omp_scaling_%j.err

module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

EXEC="./stencil_hybrid"
ARGS="-x 1000 -y 1000 -n 1000"

THREAD_COUNTS=(1 2 4 8 16 32 64 128)

echo "Threads,MinTime,MaxTime,AvgTime" > results/omp_results.csv

for t in "${THREAD_COUNTS[@]}"; do
    export OMP_NUM_THREADS=$t
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo "Running with $t threads..."
    mpirun -np 1 $EXEC $ARGS > out_omp_$t.txt
    
    min=$(grep "Min thread computation time" out_omp_$t.txt | awk '{print $5}')
    max=$(grep "Max thread computation time" out_omp_$t.txt | awk '{print $5}')
    avg=$(grep "Avg thread computation time" out_omp_$t.txt | awk '{print $5}')
    echo "$t,$min,$max,$avg" >> results/omp_results.csv
done
