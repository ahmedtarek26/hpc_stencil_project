#!binbash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=EPYC
#SBATCH --account=dssc
#SBATCH --time=000500
#SBATCH --job-name=test_stencil
#SBATCH --output=test_%j.out
#SBATCH --error=test_%j.err

# Load modules
module load gcc12.2.0
module load openmpi4.1.6--gcc--12.2.0

echo =========================================
echo ORFEO Test Run
echo =========================================
echo Node $(hostname)
echo Date $(date)
echo 

# Check if executable exists
if [ ! -f .stencil_hybrid ]; then
    echo ERROR Executable '.stencil_hybrid' not found!
    echo Please compile your code first
    echo   mpicc -fopenmp -O3 -o stencil_hybrid stencil_parallel.c -lm
    exit 1
fi

echo Executable found .stencil_hybrid
echo 

# Test with small problem
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

echo Running test with
echo   - 1 MPI task
echo   - 8 OpenMP threads
echo   - Grid 500x500
echo   - Iterations 100
echo 

srun --ntasks=1 --cpus-per-task=8 .stencil_hybrid -x 500 -y 500 -n 100

echo 
echo =========================================
echo Test completed successfully!
echo =========================================