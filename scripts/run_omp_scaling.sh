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

module purge
module load openMPI/5.0.5

EXEC="./stencil_hybrid"
ARGS="-x 2000 -y 2000 -n 1000"

# Thread counts - test up to node capacity
THREAD_COUNTS="1 2 4 8 16 32 64 128"

mkdir -p results
echo "Threads,TotalTime,CommTime,CompTime,Speedup,Efficiency" > results/omp_results.csv

BASELINE_TIME=""

for t in $THREAD_COUNTS; do
    export OMP_NUM_THREADS=$t
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    
    echo ""
    echo "====================================="
    echo "Running with $t OpenMP threads"
    echo "====================================="
    
    srun --ntasks=1 --cpus-per-task=$t $EXEC $ARGS > out_omp_${t}.txt 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Run failed for $t threads"
        continue
    fi
    
    total=$(grep "Total time" out_omp_${t}.txt | awk '{print $3}')
    comm=$(grep "Communication time" out_omp_${t}.txt | awk '{print $3}')
    comp=$(grep "Computation time" out_omp_${t}.txt | awk '{print $3}')
    
    if [ -z "$total" ]; then
        echo "WARNING: Could not extract timing for $t threads"
        total="ERROR"
        comm="ERROR"
        comp="ERROR"
        speedup="ERROR"
        efficiency="ERROR"
    else
        if [ -z "$BASELINE_TIME" ]; then
            BASELINE_TIME=$total
            speedup=1.0
            efficiency=1.0
        else
            speedup=$(echo "scale=4; $BASELINE_TIME / $total" | bc)
            efficiency=$(echo "scale=4; $speedup / $t" | bc)
        fi
    fi
    
    echo "$t,$total,$comm,$comp,$speedup,$efficiency" >> results/omp_results.csv
    echo "Completed: Total=${total}s, Speedup=${speedup}x, Efficiency=${efficiency}"
done

echo ""
echo "========================================="
echo "OpenMP scaling study complete!"
echo "Results: results/omp_results.csv"
echo "========================================="
