#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=EPYC
#SBATCH --time=00:30:00
#SBATCH --job-name=mpi_strong
#SBATCH --output=mpi_strong_%j.out
#SBATCH --error=mpi_strong_%j.err

module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

EXEC="./stencil_hybrid"
ARGS="-x 2000 -y 2000 -n 1000"

NODE_COUNTS=(1 2 4)

echo "Nodes,Tasks,Threads,CommTime,MinComp,MaxComp,AvgComp" > results/mpi_strong_results.csv

for n in "${NODE_COUNTS[@]}"; do
    export OMP_NUM_THREADS=128
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo "Running on $n nodes..."
    srun --nodes=$n --ntasks=$n --cpus-per-task=128 $EXEC $ARGS > out_mpi_strong_$n.txt
    
    comm=$(grep "Time spent in communication" out_mpi_strong_$n.txt | awk '{print $5}')
    min=$(grep "Min thread computation time" out_mpi_strong_$n.txt | awk '{print $5}')
    max=$(grep "Max thread computation time" out_mpi_strong_$n.txt | awk '{print $5}')
    avg=$(grep "Avg thread computation time" out_mpi_strong_$n.txt | awk '{print $5}')
    echo "$n,$n,128,$comm,$min,$max,$avg" >> results/mpi_strong_results.csv
done
