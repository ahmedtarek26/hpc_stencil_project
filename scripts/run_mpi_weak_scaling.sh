#!/bin/bash
#SBATCH --nodes=4              # Max allowed
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=epyc
#SBATCH --account=dssc
#SBATCH --time=01:00:00
#SBATCH --job-name=mpi_weak
#SBATCH --output=mpi_weak_%j.out
#SBATCH --error=mpi_weak_%j.err

module purge
module load openMPI/5.0.5

EXEC="./stencil_hybrid"

# Base problem size per node
BASE_SIZE=2000
ITERATIONS=1000

# Node counts - perfect squares for 2D decomposition, max 4 nodes
NODE_COUNTS="1 4"  # 1x1, 2x2 task grids

THREADS_PER_TASK=128

mkdir -p results
echo "Nodes,Tasks,ThreadsPerTask,GridX,GridY,TotalTime,CommTime,CompTime,Efficiency" > results/mpi_weak_results.csv

BASELINE_TIME=""

for n in $NODE_COUNTS; do
    # For 2D domain decomposition: scale each dimension by sqrt(n)
    SCALE=$(echo "sqrt($n)" | bc -l)
    GRID_X=$(echo "scale=0; $BASE_SIZE * $SCALE / 1" | bc)
    GRID_Y=$GRID_X
    
    export OMP_NUM_THREADS=$THREADS_PER_TASK
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo ""
    echo "====================================="
    echo "Running WEAK scaling: $n nodes ($n MPI tasks)"
    echo "Grid: ${GRID_X}x${GRID_Y}"
    echo "Iterations: $ITERATIONS"
    echo "Work per node: constant"
    echo "====================================="
    
    srun --nodes=$n --ntasks=$n --cpus-per-task=$THREADS_PER_TASK \
         $EXEC -x $GRID_X -y $GRID_Y -n $ITERATIONS > out_weak_${n}.txt 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Run failed for $n nodes"
        continue
    fi
    
    total=$(grep "Total time" out_weak_${n}.txt | awk '{print $3}')
    comm=$(grep "Communication time" out_weak_${n}.txt | awk '{print $3}')
    comp=$(grep "Computation time" out_weak_${n}.txt | awk '{print $3}')
    
    if [ -z "$total" ]; then
        echo "WARNING: Could not extract timing for $n nodes"
        continue
    fi
    
    # Calculate efficiency (ideal: time stays constant)
    if [ -z "$BASELINE_TIME" ]; then
        BASELINE_TIME=$total
        efficiency=1.0
    else
        efficiency=$(echo "scale=4; $BASELINE_TIME / $total" | bc)
    fi
    
    echo "$n,$n,$THREADS_PER_TASK,$GRID_X,$GRID_Y,$total,$comm,$comp,$efficiency" >> results/mpi_weak_results.csv
    echo "Completed: Total=${total}s, Efficiency=${efficiency}"
done

echo ""
echo "========================================="
echo "MPI WEAK scaling study complete!"
echo "Results: results/mpi_weak_results.csv"
echo "========================================="
