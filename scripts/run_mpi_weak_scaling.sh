#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=EPYC
#SBATCH --time=00:30:00
#SBATCH --job-name=mpi_weak
#SBATCH --output=mpi_weak_%j.out
#SBATCH --error=mpi_weak_%j.err

module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

EXEC="./stencil_hybrid"
BASE_SIZE=1000

NODE_COUNTS=(1 2 4)

echo "Nodes,Tasks,Threads,SizeX,SizeY,CommTime,MinComp,MaxComp,AvgComp" > results/mpi_weak_results.csv

for n in "${NODE_COUNTS[@]}"; do
    # For 2D decomposition, we scale the grid. 
    # If n=1, 1000x1000. If n=2, 1414x1414 (approx). If n=4, 2000x2000.
    # To keep it simple, let's scale one dimension for weak scaling if 1D, 
    # but the code supports 2D. Let's use 2D scaling: size = BASE_SIZE * sqrt(n)
    SIZE=$(echo "scale=0; $BASE_SIZE * sqrt($n)" | bc -l)
    
    export OMP_NUM_THREADS=128
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo "Running on $n nodes with size ${SIZE}x${SIZE}..."
    srun --nodes=$n --ntasks=$n --cpus-per-task=128 $EXEC -x $SIZE -y $SIZE -n 1000 > out_mpi_weak_$n.txt
    
    comm=$(grep "Time spent in communication" out_mpi_weak_$n.txt | awk '{print $5}')
    min=$(grep "Min thread computation time" out_mpi_weak_$n.txt | awk '{print $5}')
    max=$(grep "Max thread computation time" out_mpi_weak_$n.txt | awk '{print $5}')
    avg=$(grep "Avg thread computation time" out_mpi_weak_$n.txt | awk '{print $5}')
    echo "$n,$n,128,$SIZE,$SIZE,$comm,$min,$max,$avg" >> results/mpi_weak_results.csv
done
