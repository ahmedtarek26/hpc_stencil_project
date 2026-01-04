#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=EPYC
#SBATCH --account=dssc
#SBATCH --time=01:00:00
#SBATCH --job-name=mpi_strong
#SBATCH --output=mpi_strong_%j.out
#SBATCH --error=mpi_strong_%j.err

# Load modules
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

EXEC="./stencil_hybrid"

# Fixed problem size for strong scaling
GRID_SIZE=4000
ITERATIONS=1000

# Node counts for scaling (1, 2, 4, 8, 16)
NODE_COUNTS=(1 2 4 8 16)

# OpenMP threads per MPI task (use all cores on node)
THREADS_PER_TASK=128

# Ensure the results directory exists
mkdir -p results

echo "Nodes,Tasks,ThreadsPerTask,TotalTime,CommTime,CompTime" > results/mpi_strong_results.csv

for n in "${NODE_COUNTS[@]}"; do
    export OMP_NUM_THREADS=$THREADS_PER_TASK
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo "========================================="
    echo "Running STRONG scaling on $n nodes (${n} MPI tasks)"
    echo "Grid: ${GRID_SIZE}x${GRID_SIZE}, Iterations: ${ITERATIONS}"
    echo "========================================="
    
    # Run with n MPI tasks, one per node
    srun --nodes=$n --ntasks=$n --cpus-per-task=$THREADS_PER_TASK \
         $EXEC -x $GRID_SIZE -y $GRID_SIZE -n $ITERATIONS > out_strong_$n.txt 2>&1
    
    # Extract timing information
    total=$(grep "Total time:" out_strong_$n.txt | awk '{print $3}')
    comm=$(grep "Communication time:" out_strong_$n.txt | awk '{print $3}')
    comp=$(grep "Computation time:" out_strong_$n.txt | awk '{print $3}')
    
    # Save to CSV (calculate speedup/efficiency later in Python/Excel)
    echo "$n,$n,$THREADS_PER_TASK,$total,$comm,$comp" >> results/mpi_strong_results.csv
    
    echo "Completed: $n nodes - Total: ${total}s, Comm: ${comm}s, Comp: ${comp}s"
    echo ""
done

echo "========================================="
echo "MPI strong scaling study complete!"
echo "Results saved to: results/mpi_strong_results.csv"
echo ""
echo "To calculate speedup and efficiency:"
echo "  Speedup(n) = Time(1 node) / Time(n nodes)"
echo "  Efficiency(n) = Speedup(n) / n"
echo "========================================="