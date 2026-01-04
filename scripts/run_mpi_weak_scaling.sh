#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=EPYC
#SBATCH --account=dssc
#SBATCH --time=01:00:00
#SBATCH --job-name=mpi_weak
#SBATCH --output=mpi_weak_%j.out
#SBATCH --error=mpi_weak_%j.err

# Load modules
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

EXEC="./stencil_hybrid"

# Base problem size per node (keep workload per node constant)
BASE_SIZE=2000
ITERATIONS=1000

# Node counts for scaling (1, 2, 4, 8, 16)
NODE_COUNTS=(1 2 4 8 16)

# OpenMP threads per MPI task
THREADS_PER_TASK=128

# Ensure the results directory exists
mkdir -p results

echo "Nodes,Tasks,ThreadsPerTask,GridX,GridY,TotalTime,CommTime,CompTime,Efficiency" > results/mpi_weak_results.csv

# Store baseline time for efficiency calculation
BASELINE_TIME=""

for n in "${NODE_COUNTS[@]}"; do
    # For weak scaling: increase problem size proportionally
    # For 2D grid of tasks, scale size by sqrt(n)
    SCALE=$(echo "scale=4; sqrt($n)" | bc)
    GRID_X=$(echo "scale=0; $BASE_SIZE * $SCALE / 1" | bc)
    GRID_Y=$GRID_X
    
    export OMP_NUM_THREADS=$THREADS_PER_TASK
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo "========================================="
    echo "Running WEAK scaling on $n nodes (${n} MPI tasks)"
    echo "Grid: ${GRID_X}x${GRID_Y}, Iterations: ${ITERATIONS}"
    echo "Workload per node kept constant"
    echo "========================================="
    
    # Run with n MPI tasks, one per node
    srun --nodes=$n --ntasks=$n --cpus-per-task=$THREADS_PER_TASK \
         $EXEC -x $GRID_X -y $GRID_Y -n $ITERATIONS > out_weak_$n.txt 2>&1
    
    # Extract timing information
    total=$(grep "Total time:" out_weak_$n.txt | awk '{print $3}')
    comm=$(grep "Communication time:" out_weak_$n.txt | awk '{print $3}')
    comp=$(grep "Computation time:" out_weak_$n.txt | awk '{print $3}')
    
    # Calculate efficiency (ideal weak scaling: time stays constant)
    if [ -z "$BASELINE_TIME" ]; then
        BASELINE_TIME=$total
        efficiency=1.0
    else
        efficiency=$(echo "scale=4; $BASELINE_TIME / $total" | bc)
    fi
    
    echo "$n,$n,$THREADS_PER_TASK,$GRID_X,$GRID_Y,$total,$comm,$comp,$efficiency" >> results/mpi_weak_results.csv
    
    echo "Completed: $n nodes - Total: ${total}s, Efficiency: ${efficiency}"
    echo ""
done

echo "========================================="
echo "MPI weak scaling study complete!"
echo "Results saved to: results/mpi_weak_results.csv"
echo "========================================="