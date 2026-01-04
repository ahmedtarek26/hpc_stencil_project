# Hybrid MPI+OpenMP 2D Heat Equation Stencil Method

[![HPC](https://img.shields.io/badge/HPC-High%20Performance%20Computing-blue)](https://github.com)
[![MPI](https://img.shields.io/badge/MPI-Parallel%20Computing-green)](https://www.mpi-forum.org/)
[![OpenMP](https://img.shields.io/badge/OpenMP-Multi--threading-orange)](https://www.openmp.org/)
[![C](https://img.shields.io/badge/Language-C-lightgrey)](https://en.wikipedia.org/wiki/C_(programming_language))

A high-performance parallel implementation of the 2D heat equation using a 5-point stencil method with hybrid MPI+OpenMP parallelization. This project demonstrates efficient domain decomposition and scalability analysis on HPC systems.

## 🎯 Project Overview

This project implements a parallel solver for the 2D heat diffusion equation using:
- **MPI** for distributed-memory parallelism (domain decomposition across nodes)
- **OpenMP** for shared-memory parallelism (multi-threading within nodes)
- **5-point stencil** method for spatial discretization
- **Periodic or fixed boundary conditions**

### 📊 Key Results

| Metric | Performance |
|--------|-------------|
| **Weak Scaling Efficiency** | 110.3% (super-linear!)  |
| **Strong Scaling Efficiency (20K×20K)** | 94.5% on 4 nodes |
| **OpenMP Optimal Configuration** | 128 threads |
| **Maximum Grid Size Tested** | 20,000 × 20,000 (400M points) |

## 🚀 Features

-  Hybrid MPI+OpenMP parallelization
-  2D domain decomposition with automatic grid factorization
-  Non-blocking MPI communication (halo exchange)
-  OpenMP static scheduling for uniform workload
-  Periodic and fixed boundary conditions
-  Comprehensive performance instrumentation
-  Weak and strong scaling analysis
-  SLURM batch scripts for HPC clusters
-  Python visualization tools

## 📁 Project Structure

```
.
├── src/                              # Source code
│   ├── stencil_template_parallel.c   # Main implementation (MPI+OpenMP)
│   └── stencil_template_serial.c     
│
├── include/                          
│   ├── stencil_template_parallel.h   
│   └── stencil_template_serial.h 
│
├── scripts/                          # SLURM job submission scripts
│   ├── run_mpi_weak_scaling.sh       # Weak scaling study
│   ├── run_mpi_strong_scaling.sh     # Strong scaling study
│   ├── run_omp_scaling.sh            # OpenMP thread scaling
│   └── test_stencil.sh               # Quick test/verification
│
├── plots/                            # Visualization and analysis
│   ├── plot_mpi_weak_scaling.ipynb   # Weak scaling plots
│   ├── plot_mpi_strong_scaling.ipynb # Strong scaling plots
│   └── plot_omp_scaling.ipynb        # OpenMP scaling plots
│
├── results/                          # Output data (CSV files)
│   ├── mpi_weak_results.csv          # Weak scaling timings
│   ├── mpi_strong_results.csv        # Strong scaling timings
│   └── omp_results.csv               # OpenMP timings
│
├── README.md                         # This file
└── LICENSE                           # Project license
```

## 🛠️ Requirements

### Software Dependencies
- **Compiler**: GCC 12.2.0+ (with OpenMP support)
- **MPI Implementation**: OpenMPI 4.1.6+ or MPICH
- **SLURM**: For HPC job scheduling (optional)
- **Python 3.8+**: For visualization (optional)

### Python Libraries (for plotting)
```bash
pip install numpy matplotlib pandas seaborn jupyter
```

## 🔨 Compilation

### Basic Compilation
```bash
cd src/
mpicc -fopenmp -O3 -o stencil_hybrid stencil_template_parallel.c -lm
```

### With Debugging
```bash
mpicc -fopenmp -g -O0 -o stencil_hybrid_debug stencil_template_parallel.c -lm
```

### Architecture-Specific Optimization
```bash
# For Intel CPUs
mpicc -fopenmp -O3 -march=native -mtune=native -o stencil_hybrid stencil_template_parallel.c -lm

# For AMD EPYC CPUs
mpicc -fopenmp -O3 -march=znver3 -mtune=znver3 -o stencil_hybrid stencil_template_parallel.c -lm
```

## 🏃 Usage

### Command-Line Options

```bash
./stencil_hybrid [options]

Options:
  -x <size>    Grid width (default: 1000)
  -y <size>    Grid height (default: 1000)
  -n <iters>   Number of iterations (default: 100)
  -e <sources> Number of energy sources (default: 1)
  -E <energy>  Energy per source (default: 1.0)
  -p <0|1>     Periodic boundaries (default: 0)
  -o <0|1>     Print energy at every step (default: 0)
  -v <level>   Verbosity level (default: 0)
  -h           Display help
```

### Running Locally

#### Single Node, Multiple Threads
```bash
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

mpirun -np 1 ./stencil_hybrid -x 2000 -y 2000 -n 1000
```

#### Multiple MPI Tasks
```bash
export OMP_NUM_THREADS=16

mpirun -np 4 ./stencil_hybrid -x 4000 -y 4000 -n 1000
```

### Running on HPC Cluster (SLURM)

#### Test Run
```bash
cd scripts/
sbatch test_stencil.sh
```

#### Weak Scaling Study
```bash
sbatch run_mpi_weak_scaling.sh
```

#### Strong Scaling Study
```bash
sbatch run_mpi_strong_scaling.sh
```

#### OpenMP Scaling Study
```bash
sbatch run_omp_scaling.sh
```

## 📈 Performance Results

### Weak Scaling (Constant Work per Node)
- **Problem size per node**: ~2000×2000
- **Nodes tested**: 1, 2, 4
- **Efficiency**: 105.2% → **110.3%** (super-linear scaling!)
- **Reason**: Improved cache behavior with larger total problem size

### Strong Scaling (Fixed Problem Size)
#### Grid 10,000×10,000
- **Nodes**: 1, 2, 4
- **Speedup**: 1.00x → 1.40x → 1.90x
- **Efficiency**: 100% → 70.0% → 47.5%
- **Communication overhead**: 6.4% → 13.5% → 21.0%

#### Grid 20,000×20,000
- **Nodes**: 1, 2, 4
- **Speedup**: 1.00x → 1.96x → 3.78x
- **Efficiency**: 100% → 97.8% → **94.5%**
- **Communication overhead**: 0.2% → 1.7% → 5.0%

### OpenMP Scaling (Single Node)
- **Optimal threads**: 128
- **Minimum time**: 0.139s at 128 threads
- **Key finding**: Static scheduling outperforms dynamic/guided

## 🔬 Algorithm Details

### 5-Point Stencil Formula
```
u[i,j]^(n+1) = (u[i-1,j]^n + u[i+1,j]^n + u[i,j-1]^n + u[i,j+1]^n) / 4
```

### Domain Decomposition
- **2D Cartesian grid** of MPI processes
- **Halo exchange** using non-blocking MPI_Isend/MPI_Irecv
- **Automatic factorization** to create square-like process grid

### Communication Pattern
```
       [North Neighbor]
              |
[West] -- [Current] -- [East]
              |
       [South Neighbor]
```

## 📊 Visualization

Generate plots using the provided Jupyter notebooks:

```bash
cd plots/
jupyter notebook plot_mpi_weak_scaling.ipynb
jupyter notebook plot_mpi_strong_scaling.ipynb  
jupyter notebook plot_omp_scaling.ipynb
```

## 🎓 Academic Context

This project was developed as part of the **Foundations of High Performance Computing** course at **Università di Trieste** (2024-2025), taught by Prof. Luca Tornatore.

### Platform
- **System**: Leonardo/ORFEO Supercomputer
- **Nodes**: Dual-socket AMD EPYC with 128 cores/node
- **Interconnect**: High-speed network for MPI communication

## 🤝 Contributing

Contributions are welcome! Areas for improvement:
- Implement communication-computation overlap
- Add GPU acceleration support
- Implement additional boundary conditions
- Add checkpoint/restart capability
- Optimize for different network topologies



## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.


## 📧 Contact

For questions or collaboration:
- **Author**: Ahmed Tarek
- **Email**: ahmedtarek2632@gmail.com

---

**⭐ If you find this project useful, please consider giving it a star!**
