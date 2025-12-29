# High-Performance Computing Project: 2D Stencil Method on ORFEO

**Author:** Manus AI
**Date:** December 28, 2025

## 1. Introduction and Project Goal

This document outlines the implementation and execution plan for the High-Performance Computing assignment, which focuses on solving the 2D Heat Equation using a **hybrid MPI+OpenMP stencil method** [1]. The primary goal is to implement an efficient parallel solution and conduct a thorough performance analysis, specifically focusing on **OpenMP scaling**, **MPI strong scaling**, and **MPI weak scaling** on the ORFEO cluster [2].

The provided solution is based on the course template and incorporates the parallelization strategies successfully implemented by a peer [3].

## 2. Implementation Details

The core of the project is the `stencil_hybrid` executable, which uses a hybrid parallel programming model:

*   **MPI (Message Passing Interface):** Used for domain decomposition, dividing the 2D grid into smaller patches, with each patch assigned to a separate MPI process. This handles inter-node and inter-process communication.
*   **OpenMP (Open Multi-Processing):** Used for shared-memory parallelization within each MPI process, distributing the computation of the local patch across multiple threads.

### 2.1. Parallelization Strategy

The implementation uses a standard 5-point stencil for the heat diffusion update.

| Component | Parallelization Technique | Key Feature |
| :--- | :--- | :--- |
| **Domain Decomposition** | MPI | 2D decomposition is used to minimize the surface-to-volume ratio, thus reducing communication overhead. |
| **Halo Exchange** | MPI Non-Blocking Communication (`MPI_Isend`, `MPI_Irecv`) | Non-blocking calls are used to overlap communication with computation. The `MPI_Waitall` is called after the local computation to ensure all halo data is received before the next iteration. |
| **Local Computation** | OpenMP | The main stencil loop (`update_plane` function) is parallelized using `#pragma omp parallel for collapse(2)`. The `collapse(2)` clause parallelizes the nested loops over the local x and y dimensions, improving load balance and reducing OpenMP overhead. |
| **Timing** | MPI and OpenMP Timers | The code includes instrumentation to measure total wall-clock time, MPI communication time, and per-thread OpenMP computation time for detailed analysis. |

## 3. Preparing the Project for ORFEO

The project is packaged in a zip file named `hpc_stencil_project_orfeo.zip`.

### 3.1. Required Modules

The ORFEO cluster uses the **LMOD** system for managing software environments. Based on the cluster documentation and peer's successful runs, the following modules are required for compilation and execution:

```bash
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0
```

### 3.2. Compilation

The provided `Makefile` automates the compilation process.

1.  **Transfer the project:** Upload the `hpc_stencil_project_orfeo.zip` file to your home directory on ORFEO and unzip it.
    ```bash
    unzip hpc_stencil_project_orfeo.zip
    cd hpc_stencil_project
    ```
2.  **Load modules and compile:**
    ```bash
    module load gcc/12.2.0
    module load openmpi/4.1.6--gcc--12.2.0
    make
    ```
    This will create the executable file `./stencil_hybrid`.

## 4. Scaling Experiments on ORFEO (SLURM Scripts)

Three SLURM batch scripts are provided in the `scripts/` directory to automate the required scaling experiments on the **EPYC** partition of the ORFEO cluster. The scripts are configured to use the **AMD EPYC 7H12** nodes, which have 128 cores per node.

### 4.1. OpenMP Scaling (`scripts/run_omp_scaling.sh`)

This script tests the performance of the OpenMP parallelization on a single node with a fixed problem size (1000x1000 grid). It varies the number of threads from 1 up to 128.

| SLURM Parameter | Value | Description |
| :--- | :--- | :--- |
| `--nodes` | 1 | Single node test. |
| `--ntasks-per-node` | 1 | One MPI process per node (the master process). |
| `--cpus-per-task` | 128 | Maximum number of cores available for OpenMP threads. |
| `OMP_NUM_THREADS` | 1, 2, 4, ..., 128 | The number of OpenMP threads is varied in a loop. |

**Execution:** `sbatch scripts/run_omp_scaling.sh`

### 4.2. MPI Strong Scaling (`scripts/run_mpi_strong_scaling.sh`)

This script tests how the runtime decreases as the number of nodes (and thus total cores) increases, while keeping the **total problem size fixed** (2000x2000 grid). It uses 128 threads per MPI task, effectively running one MPI task per node.

| SLURM Parameter | Value | Description |
| :--- | :--- | :--- |
| `--nodes` | 1, 2, 4 | Number of nodes is varied in a loop. |
| `--ntasks` | 1, 2, 4 | Number of MPI tasks is equal to the number of nodes. |
| `OMP_NUM_THREADS` | 128 | All cores on the node are used by OpenMP threads. |

**Execution:** `sbatch scripts/run_mpi_strong_scaling.sh`

### 4.3. MPI Weak Scaling (`scripts/run_mpi_weak_scaling.sh`)

This script tests how the runtime changes as both the number of nodes and the **problem size per node** increase, keeping the workload per core constant. The total problem size scales with the number of nodes (e.g., 1000x1000 for 1 node, 2000x2000 for 4 nodes).

| SLURM Parameter | Value | Description |
| :--- | :--- | :--- |
| `--nodes` | 1, 2, 4 | Number of nodes is varied in a loop. |
| `--ntasks` | 1, 2, 4 | Number of MPI tasks is equal to the number of nodes. |
| `Problem Size` | Scales with $\sqrt{N_{nodes}}$ | Workload per node remains constant. |

**Execution:** `sbatch scripts/run_mpi_weak_scaling.sh`

## 5. Analyzing Results

Each script will generate a CSV file containing the performance metrics:

*   `omp_results.csv`
*   `mpi_strong_results.csv`
*   `mpi_weak_results.csv`

These files contain the necessary timing data (communication time, min/max/avg computation time) to calculate **Speedup** and **Efficiency** for your final report, as required by the assignment [1].

---

## References

[1] [slides\_assignment\_2025.v1.2\_July\_26th.pdf](file:///home/ubuntu/upload/slides_assignment_2025.v1.2_July_26th.pdf)
[2] [ORFEO documentation](https://orfeo-doc.areasciencepark.it/)
[3] [zahra-ynp/HPCProject-StencilMethod](https://github.com/zahra-ynp/HPCProject-StencilMethod/tree/main)
