#!/bin/bash
#SBATCH --nodes=4              # Max allowed for dssc account
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=epyc       # lowercase
#SBATCH --account=dssc
#SBATCH --time=01:00:00
#SBATCH --job-name=mpi_strong
#SBATCH --output=mpi_strong_%j.out
#SBATCH --error=mpi_strong_%j.err

module purge
module load openMPI/5.0.5

EXEC="./stencil_hybrid"

# Fixed problem size for strong scaling
GRIDSIZE=4000
ITERATIONS=1000

# Node counts - LIMITED TO 4 NODES MAX on EPYC/dssc
NODE_COUNTS="1 2 4"

# OpenMP threads per MPI task (adjust based on actual cores per EPYC node)
THREADS_PER_TASK=128

mkdir -p results
echo "Nodes,Tasks,ThreadsPerTask,TotalTime,CommTime,CompTime,Speedup,Efficiency" > results/mpi_strong_results.csv

BASELINE_TIME=""

for n in $NODE_COUNTS; do
    export OMP_NUM_THREADS=$THREADS_PER_TASK
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo ""
    echo "====================================="
    echo "Running STRONG scaling: $n nodes ($n MPI tasks)"
    echo "Grid: ${GRIDSIZE}x${GRIDSIZE}"
    echo "Iterations: $ITERATIONS"
    echo "Threads per task: $THREADS_PER_TASK"
    echo "====================================="
    
    srun --nodes=$n --ntasks=$n --cpus-per-task=$THREADS_PER_TASK \
         $EXEC -x $GRIDSIZE -y $GRIDSIZE -n $ITERATIONS > out_strong_${n}.txt 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Run failed for $n nodes"
        continue
    fi
    
    total=$(grep "Total time" out_strong_${n}.txt | awk '{print $3}')
    comm=$(grep "Communication time" out_strong_${n}.txt | awk '{print $3}')
    comp=$(grep "Computation time" out_strong_${n}.txt | awk '{print $3}')
    
    if [ -z "$total" ]; then
        echo "WARNING: Could not extract timing for $n nodes"
        continue
    fi
    
    # Store baseline for speedup calculation
    if [ -z "$BASELINE_TIME" ]; then
        BASELINE_TIME=$total
        speedup=1.0
        efficiency=1.0
    else
        speedup=$(echo "scale=4; $BASELINE_TIME / $total" | bc)
        efficiency=$(echo "scale=4; $speedup / $n" | bc)
    fi
    
    echo "$n,$n,$THREADS_PER_TASK,$total,$comm,$comp,$speedup,$efficiency" >> results/mpi_strong_results.csv
    echo "Completed: Total=${total}s, Speedup=${speedup}x, Efficiency=${efficiency}"
done

echo ""
echo "========================================="
echo "MPI STRONG scaling study complete!"
echo "Results: results/mpi_strong_results.csv"
echo "========================================="
