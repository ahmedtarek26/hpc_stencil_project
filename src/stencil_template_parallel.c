/*
 * Fixed and Complete Parallel Stencil Implementation
 * MPI + OpenMP Hybrid Parallelization
 * 5-point stencil for heat equation
 */

#include "stencil_template_parallel.h"

// --- GLOBAL VARIABLE DEFINITIONS ---
int g_n_omp_threads = 1;
double* g_per_thread_comp_time = NULL;
__thread double thread_local_comp_time = 0.0;

// ------------------------------------------------------------------
// FUNCTION PROTOTYPES (if not in header)
// ------------------------------------------------------------------
inline int inject_energy(const int periodic, const int Nsources, const vec2_t *Sources,
                         const double energy, plane_t *plane, const vec2_t N);
inline int update_plane(const int periodic, const vec2_t N, const plane_t *oldplane, plane_t *newplane);
inline int get_total_energy(plane_t *plane, double *energy);
int memory_release(plane_t *planes, buffers_t *buffers, vec2_t *sources_local);
int memory_allocate(const uint neighbours[4], const vec2_t N, buffers_t *buffers_ptr, plane_t *planes_ptr);
int output_energy_stat(int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm);
int initialize(MPI_Comm *Comm, int Me, int Ntasks, int argc, char **argv, vec2_t *S, vec2_t *N,
               int *periodic, int *output_energy_stat_perstep, uint *neighbours, int *Niterations,
               int *Nsources, int *Nsources_local, vec2_t **Sources_local, double *energy_per_source,
               plane_t *planes, buffers_t *buffers);
uint simple_factorization(uint A, int *Nfactors, uint **factors);
int initialize_sources(int Me, int Ntasks, MPI_Comm *Comm, vec2_t mysize, int Nsources,
                       int *Nsources_local, vec2_t **Sources);

// ------------------------------------------------------------------
// MAIN
// ------------------------------------------------------------------
int main(int argc, char **argv)
{
  MPI_Comm myCOMM_WORLD;
  int  Rank, Ntasks;
  uint neighbours[4];

  int  Niterations;
  int  periodic;
  vec2_t S, N;
  
  int      Nsources;
  int      Nsources_local;
  vec2_t  *Sources_local;
  double   energy_per_source;

  plane_t   planes[2];  
  buffers_t buffers[2];
  
  int output_energy_stat_perstep;
  
  /* Initialize MPI environment */
  {
    int level_obtained;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &level_obtained);
    if (level_obtained < MPI_THREAD_FUNNELED) {
      printf("MPI_thread level obtained is %d instead of %d\n",
             level_obtained, MPI_THREAD_FUNNELED);
      MPI_Finalize();
      exit(1);
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
    MPI_Comm_dup(MPI_COMM_WORLD, &myCOMM_WORLD);
  }
  
  /* Setup for per-thread timing */
  #pragma omp parallel
  {
    #pragma omp master
    { 
      g_n_omp_threads = omp_get_num_threads(); 
    }
  }
  g_per_thread_comp_time = (double*)calloc(g_n_omp_threads, sizeof(double));

  /* Argument checking and setting */
  int ret = initialize(&myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, 
                       &output_energy_stat_perstep, neighbours, &Niterations,
                       &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
                       &planes[0], &buffers[0]);

  if (ret) {
    printf("Task %d is opting out with termination code %d\n", Rank, ret);
    MPI_Finalize();
    return 0;
  }
  
  int current = OLD;
  double t_start, t_end;
  double total_comm_time = 0.0;
  double total_comp_time = 0.0;
  
  t_start = MPI_Wtime();

  /* Main iteration loop */
  for (int iter = 0; iter < Niterations; ++iter) {
    double section_start_time;
 
    /* Inject new energy from sources */
    inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
    
    section_start_time = MPI_Wtime();

    /* ====================================== */
    /* COMMUNICATION PHASE: Halo Exchange     */
    /* ====================================== */

    double* current_plane = planes[current].data;
    const int sizex = planes[current].size[_x_];
    const int sizey = planes[current].size[_y_];
    const int full_sizex = sizex + 2;
    const int full_sizey = sizey + 2;  // FIX: Define full_sizey

    /* Fill the buffers and/or make the buffers' pointers point to correct positions */
    
    // North and South buffers can point directly to the data (contiguous lines)
    buffers[SEND][NORTH] = current_plane + full_sizex;  // First inner row
    buffers[SEND][SOUTH] = current_plane + sizey * full_sizex;  // Last inner row

    buffers[RECV][NORTH] = current_plane;  // Halo above first row
    buffers[RECV][SOUTH] = current_plane + (sizey + 1) * full_sizex;  // Halo below last row

    // East and West: Fill send buffers (non-contiguous data)
    #pragma omp parallel for
    for (int j = 1; j <= sizey; j++) {
      buffers[SEND][WEST][j-1] = current_plane[j * full_sizex + 1];  // First inner column
      buffers[SEND][EAST][j-1] = current_plane[j * full_sizex + sizex];  // Last inner column
    }

    /* Perform the halo communications using non-blocking Isend/Irecv */
    MPI_Request reqs[8];
    int req_idx = 0;

    // Send/Recv North
    if (neighbours[NORTH] != MPI_PROC_NULL) {
      MPI_Isend(buffers[SEND][NORTH], sizex, MPI_DOUBLE, neighbours[NORTH], 0, myCOMM_WORLD, &reqs[req_idx++]);
      MPI_Irecv(buffers[RECV][NORTH], sizex, MPI_DOUBLE, neighbours[NORTH], 1, myCOMM_WORLD, &reqs[req_idx++]);
    }

    // Send/Recv South
    if (neighbours[SOUTH] != MPI_PROC_NULL) {
      MPI_Isend(buffers[SEND][SOUTH], sizex, MPI_DOUBLE, neighbours[SOUTH], 1, myCOMM_WORLD, &reqs[req_idx++]);
      MPI_Irecv(buffers[RECV][SOUTH], sizex, MPI_DOUBLE, neighbours[SOUTH], 0, myCOMM_WORLD, &reqs[req_idx++]);
    }

    // Send/Recv West
    if (neighbours[WEST] != MPI_PROC_NULL) {
      MPI_Isend(buffers[SEND][WEST], sizey, MPI_DOUBLE, neighbours[WEST], 2, myCOMM_WORLD, &reqs[req_idx++]);
      MPI_Irecv(buffers[RECV][WEST], sizey, MPI_DOUBLE, neighbours[WEST], 3, myCOMM_WORLD, &reqs[req_idx++]);
    }

    // Send/Recv East
    if (neighbours[EAST] != MPI_PROC_NULL) {
      MPI_Isend(buffers[SEND][EAST], sizey, MPI_DOUBLE, neighbours[EAST], 3, myCOMM_WORLD, &reqs[req_idx++]);
      MPI_Irecv(buffers[RECV][EAST], sizey, MPI_DOUBLE, neighbours[EAST], 2, myCOMM_WORLD, &reqs[req_idx++]);
    }

    // Wait for all communications to complete
    MPI_Waitall(req_idx, reqs, MPI_STATUSES_IGNORE);

    total_comm_time += MPI_Wtime() - section_start_time;

    /* Copy the received haloes data for East/West */
    #pragma omp parallel for
    for (int j = 1; j <= sizey; j++) {
      current_plane[j * full_sizex] = buffers[RECV][WEST][j-1];  // Left halo
      current_plane[j * full_sizex + sizex + 1] = buffers[RECV][EAST][j-1];  // Right halo
    }

    /* If non-periodic and at boundary, set halo to 0 (infinite sink) */
    if (!periodic) {
      if (neighbours[NORTH] == MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = 0; i < full_sizex; i++) 
          current_plane[i] = 0.0;
      }
      if (neighbours[SOUTH] == MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = 0; i < full_sizex; i++) 
          current_plane[(sizey + 1) * full_sizex + i] = 0.0;
      }
      if (neighbours[WEST] == MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = 0; j < full_sizey; j++) 
          current_plane[j * full_sizex] = 0.0;
      }
      if (neighbours[EAST] == MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = 0; j < full_sizey; j++) 
          current_plane[j * full_sizex + sizex + 1] = 0.0;
      }
    }

    /* ====================================== */
    /* COMPUTATION PHASE: Update grid points */
    /* ====================================== */

    section_start_time = MPI_Wtime();
    update_plane(periodic, N, &planes[current], &planes[!current]);
    total_comp_time += MPI_Wtime() - section_start_time;

    /* Output if needed */
    if (output_energy_stat_perstep)
      output_energy_stat(iter, &planes[!current], (iter+1) * Nsources * energy_per_source, Rank, &myCOMM_WORLD);
      
    /* Swap plane indexes for the new iteration */
    current = !current;
  }
  
  t_end = MPI_Wtime() - t_start;

  /* Aggregate and print timings (reduce across ranks) */
  double max_comm_time, max_comp_time, max_total_time;
  MPI_Reduce(&total_comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
  MPI_Reduce(&total_comp_time, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
  MPI_Reduce(&t_end, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);

  if (Rank == 0) {
    printf("========================================\n");
    printf("Performance Results:\n");
    printf("========================================\n");
    printf("Total time:         %f seconds\n", max_total_time);
    printf("Communication time: %f seconds (%.2f%%)\n", max_comm_time, 100.0 * max_comm_time / max_total_time);
    printf("Computation time:   %f seconds (%.2f%%)\n", max_comp_time, 100.0 * max_comp_time / max_total_time);
    printf("========================================\n");
  }

  /* Output final energy */
  output_energy_stat(-1, &planes[current], Niterations * Nsources * energy_per_source, Rank, &myCOMM_WORLD);

  /* Cleanup */
  memory_release(planes, buffers, Sources_local);
  free(g_per_thread_comp_time);
  MPI_Finalize();
  return 0;
}

// ------------------------------------------------------------------
// INJECT ENERGY
// ------------------------------------------------------------------
inline int inject_energy(const int periodic, const int Nsources, const vec2_t *Sources,
                         const double energy, plane_t *plane, const vec2_t N)
{
  const uint sizex = plane->size[_x_] + 2;
  double * restrict data = plane->data;
    
  #define IDX(i, j) ((j) * sizex + (i))
  
  for (int s = 0; s < Nsources; s++) {
    int x = Sources[s][_x_];
    int y = Sources[s][_y_];
    
    data[IDX(x, y)] += energy;
    
    if (periodic) {
      if (N[_x_] == 1) {
        // Propagate x boundaries if single column of tasks
        if (x == 1) 
          data[IDX(plane->size[_x_] + 1, y)] += energy;
        if (x == plane->size[_x_]) 
          data[IDX(0, y)] += energy;
      }
      
      if (N[_y_] == 1) {
        // Propagate y boundaries if single row of tasks
        if (y == 1) 
          data[IDX(x, plane->size[_y_] + 1)] += energy;
        if (y == plane->size[_y_]) 
          data[IDX(x, 0)] += energy;
      }
    }                
  }
  
  #undef IDX
  return 0;
}

// ------------------------------------------------------------------
// UPDATE PLANE - 5-point stencil with OpenMP
// ------------------------------------------------------------------
inline int update_plane(const int periodic, const vec2_t N, const plane_t *oldplane, plane_t *newplane)
{
  uint fxsize = oldplane->size[_x_] + 2;
  uint fysize = oldplane->size[_y_] + 2;
  
  uint xsize = oldplane->size[_x_];
  uint ysize = oldplane->size[_y_];
  
  #define IDX(i, j) ((j) * fxsize + (i))
  
  double * restrict old = oldplane->data;
  double * restrict new = newplane->data;
  
  /* 5-point stencil update with OpenMP parallelization */
  #pragma omp parallel for collapse(2)
  for (uint j = 1; j <= ysize; j++) {
    for (uint i = 1; i <= xsize; i++) {
      new[IDX(i, j)] = old[IDX(i, j)] / 2.0 + 
                       (old[IDX(i-1, j)] + old[IDX(i+1, j)] +
                        old[IDX(i, j-1)] + old[IDX(i, j+1)]) / 4.0 / 2.0;
    }
  }

  /* Handle periodic boundaries for single task dimensions */
  if (periodic) {
    if (N[_x_] == 1) {
      #pragma omp parallel for
      for (uint j = 1; j <= ysize; j++) {
        new[IDX(0, j)] = new[IDX(xsize, j)];
        new[IDX(xsize + 1, j)] = new[IDX(1, j)];
      }
    }
  
    if (N[_y_] == 1) {
      #pragma omp parallel for
      for (uint i = 1; i <= xsize; i++) {
        new[IDX(i, 0)] = new[IDX(i, ysize)];
        new[IDX(i, ysize + 1)] = new[IDX(i, 1)];
      }
    }
  }

  #undef IDX
  return 0;
}

// ------------------------------------------------------------------
// GET TOTAL ENERGY
// ------------------------------------------------------------------
inline int get_total_energy(plane_t *plane, double *energy)
{
  const int xsize = plane->size[_x_];
  const int ysize = plane->size[_y_];
  const int fsize = xsize + 2;

  double * restrict data = plane->data;
  
  #define IDX(i, j) ((j) * fsize + (i))

  #if defined(LONG_ACCURACY)    
  long double totenergy = 0;
  #else
  double totenergy = 0;    
  #endif

  #pragma omp parallel for collapse(2) reduction(+:totenergy)
  for (int j = 1; j <= ysize; j++) {
    for (int i = 1; i <= xsize; i++) {
      totenergy += data[IDX(i, j)];
    }
  }
  
  #undef IDX

  *energy = (double)totenergy;
  return 0;
}

// ------------------------------------------------------------------
// MEMORY RELEASE
// ------------------------------------------------------------------
int memory_release(plane_t *planes, buffers_t *buffers, vec2_t *sources_local)
{
  if (planes != NULL) {
    if (planes[OLD].data != NULL)
      free(planes[OLD].data);
    
    if (planes[NEW].data != NULL)
      free(planes[NEW].data);
  }

  // Free buffers
  if (buffers != NULL) {
    if (buffers[SEND][EAST] != NULL) free(buffers[SEND][EAST]);
    if (buffers[RECV][EAST] != NULL) free(buffers[RECV][EAST]);
    if (buffers[SEND][WEST] != NULL) free(buffers[SEND][WEST]);
    if (buffers[RECV][WEST] != NULL) free(buffers[RECV][WEST]);
    // North/South are pointers into the plane data, no free needed
  }

  if (sources_local != NULL) 
    free(sources_local);
      
  return 0;
}

// ------------------------------------------------------------------
// MEMORY ALLOCATE
// ------------------------------------------------------------------
int memory_allocate(const uint neighbours[4], const vec2_t N, buffers_t *buffers_ptr, plane_t *planes_ptr)
{
  if (planes_ptr == NULL) {
    fprintf(stderr, "Error: Invalid planes_ptr.\n");
    return 1;
  }

  if (buffers_ptr == NULL) {
    fprintf(stderr, "Error: Invalid buffers_ptr.\n");
    return 1;
  }
    
  /* Allocate memory for planes with halo (frame) */
  unsigned int frame_size = (planes_ptr[OLD].size[_x_] + 2) * (planes_ptr[OLD].size[_y_] + 2);

  planes_ptr[OLD].data = (double*)malloc(frame_size * sizeof(double));
  if (planes_ptr[OLD].data == NULL) {
    perror("malloc failed for OLD plane");
    return 1;
  }
  memset(planes_ptr[OLD].data, 0, frame_size * sizeof(double));

  planes_ptr[NEW].data = (double*)malloc(frame_size * sizeof(double));
  if (planes_ptr[NEW].data == NULL) {
    perror("malloc failed for NEW plane");
    return 1;
  }
  memset(planes_ptr[NEW].data, 0, frame_size * sizeof(double));

  /* Allocate buffers for East/West (sizey doubles) */
  int sizey = planes_ptr[OLD].size[_y_];

  (*buffers_ptr)[SEND][EAST] = (double*)malloc(sizey * sizeof(double));
  (*buffers_ptr)[RECV][EAST] = (double*)malloc(sizey * sizeof(double));
  (*buffers_ptr)[SEND][WEST] = (double*)malloc(sizey * sizeof(double));
  (*buffers_ptr)[RECV][WEST] = (double*)malloc(sizey * sizeof(double));

  if (!(*buffers_ptr)[SEND][EAST] || !(*buffers_ptr)[RECV][EAST] ||
      !(*buffers_ptr)[SEND][WEST] || !(*buffers_ptr)[RECV][WEST]) {
    perror("malloc failed for buffers");
    return 1;
  }

  /* North/South buffers will point to plane data (set dynamically in main loop) */
  (*buffers_ptr)[SEND][NORTH] = NULL;
  (*buffers_ptr)[RECV][NORTH] = NULL;
  (*buffers_ptr)[SEND][SOUTH] = NULL;
  (*buffers_ptr)[RECV][SOUTH] = NULL;

  return 0;
}

// ------------------------------------------------------------------
// OUTPUT ENERGY STAT
// ------------------------------------------------------------------
int output_energy_stat(int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm)
{
  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy(plane, &system_energy);
  
  MPI_Reduce(&system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm);
  
  if (Me == 0) {
    if (step >= 0)
      printf("[ step %4d ] ", step);
    
    printf("total injected energy is %g, "
           "system energy is %g "
           "(in avg %g per grid point)\n",
           budget, tot_system_energy,
           tot_system_energy / (plane->size[_x_] * plane->size[_y_]));
  }
  
  return 0;
}

// ------------------------------------------------------------------
// INITIALIZE
// ------------------------------------------------------------------
int initialize(MPI_Comm *Comm, int Me, int Ntasks, int argc, char **argv,
               vec2_t *S, vec2_t *N, int *periodic, int *output_energy_stat_perstep,
               uint *neighbours, int *Niterations, int *Nsources, int *Nsources_local,
               vec2_t **Sources_local, double *energy_per_source,
               plane_t *planes, buffers_t *buffers)
{
  int ret = 0;
  int verbose = 0;

  vec2_t Grid;

  /* Argument parsing */
  int opt;

  /* Default values */
  (*S)[_x_] = 1000;
  (*S)[_y_] = 1000;
  *Nsources = 1;
  *energy_per_source = 1.0;
  *Niterations = 100;
  *periodic = 0;
  *output_energy_stat_perstep = 0;

  while ((opt = getopt(argc, argv, "x:y:e:E:n:p:o:v:h")) != -1) {
    switch (opt) {
      case 'x':
        (*S)[_x_] = atoi(optarg);
        break;
      case 'y':
        (*S)[_y_] = atoi(optarg);
        break;
      case 'e':
        *Nsources = atoi(optarg);
        break;
      case 'E':
        *energy_per_source = atof(optarg);
        break;
      case 'n':
        *Niterations = atoi(optarg);
        break;
      case 'p':
        *periodic = atoi(optarg);
        break;
      case 'o':
        *output_energy_stat_perstep = atoi(optarg);
        break;
      case 'v':
        verbose = atoi(optarg);
        break;
      case 'h':
        if (Me == 0) {
          printf("usage: %s [options]\n"
                 "options (overriding the default values):\n"
                 "  -x    x size of the plate [1000]\n"
                 "  -y    y size of the plate [1000]\n"
                 "  -e    how many energy sources on the plate [1]\n"
                 "  -E    energy per source [1.0]\n"
                 "  -n    how many iterations [100]\n"
                 "  -p    whether periodic boundaries apply [0 = false]\n"
                 "  -o    whether to print energy budget at every step [0 = false]\n"
                 "  -v    verbosity level [0]\n",
                 argv[0]);
        }
        MPI_Finalize();
        exit(0);
        break;
      case ':':
        if (Me == 0) printf("option -%c requires an argument\n", optopt);
        MPI_Finalize();
        exit(1);
        break;
      case '?':
        if (Me == 0) printf("-------- help unavailable ----------\n");
        MPI_Finalize();
        exit(1);
        break;
    }
  }

  /* Parameter validation */
  if ((*S)[_x_] <= 0 || (*S)[_y_] <= 0 || *Nsources <= 0 || *Niterations <= 0) {
    if (Me == 0) printf("Error: Invalid parameters.\n");
    MPI_Finalize();
    exit(1);
  }

  /* Domain decomposition - find a good 2D grid factorization */
  int Nfactors;
  uint *factors;
  simple_factorization(Ntasks, &Nfactors, &factors);

  // Try to make the grid as square as possible
  Grid[_x_] = 1;
  Grid[_y_] = Ntasks;  // Default 1D
  
  int sqrt_n = (int)sqrt(Ntasks);
  for (int i = sqrt_n; i > 0; i--) {
    if (Ntasks % i == 0) {
      Grid[_x_] = i;
      Grid[_y_] = Ntasks / i;
      break;
    }
  }

  free(factors);

  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];

  /* Compute coordinates of this task in the grid */
  int X = Me % Grid[_x_];
  int Y = Me / Grid[_x_];

  /* Set neighbors (considering periodic boundaries if requested) */
  neighbours[NORTH] = (Y > 0) ? Me - Grid[_x_] : 
                      (*periodic ? Me + (Grid[_y_] - 1) * Grid[_x_] : MPI_PROC_NULL);
  neighbours[SOUTH] = (Y < Grid[_y_] - 1) ? Me + Grid[_x_] : 
                      (*periodic ? Me - (Grid[_y_] - 1) * Grid[_x_] : MPI_PROC_NULL);
  neighbours[WEST] = (X > 0) ? Me - 1 : 
                     (*periodic ? Me + (Grid[_x_] - 1) : MPI_PROC_NULL);
  neighbours[EAST] = (X < Grid[_x_] - 1) ? Me + 1 : 
                     (*periodic ? Me - (Grid[_x_] - 1) : MPI_PROC_NULL);

  /* Compute local size (distribute grid points as evenly as possible) */
  vec2_t mysize;
  uint s = (*S)[_x_] / Grid[_x_];
  uint r = (*S)[_x_] % Grid[_x_];
  mysize[_x_] = s + (X < r ? 1 : 0);
  
  s = (*S)[_y_] / Grid[_y_];
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + (Y < r ? 1 : 0);

  planes[OLD].size[_x_] = mysize[_x_];
  planes[OLD].size[_y_] = mysize[_y_];
  planes[NEW].size[_x_] = mysize[_x_];
  planes[NEW].size[_y_] = mysize[_y_];

  if (verbose > 0 && Me == 0) {
    printf("Grid decomposition: %d x %d MPI tasks\n", Grid[_x_], Grid[_y_]);
    printf("Global grid size: %d x %d\n", (*S)[_x_], (*S)[_y_]);
  }

  /* Allocate memory */
  ret = memory_allocate(neighbours, *N, buffers, planes);
  if (ret) return ret;

  /* Initialize sources */
  ret = initialize_sources(Me, Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local);
  if (ret) return ret;

  return 0;
}

// ------------------------------------------------------------------
// SIMPLE FACTORIZATION
// ------------------------------------------------------------------
uint simple_factorization(uint A, int *Nfactors, uint **factors)
{
  int N = 0;
  int f = 2;
  uint _A_ = A;

  /* Count factors */
  while (f * f <= _A_) {
    while (_A_ % f == 0) {
      N++;
      _A_ /= f;
    }
    f++;
  }
  if (_A_ > 1) N++;  // Last prime

  *Nfactors = N;
  uint *_factors_ = (uint*)malloc(N * sizeof(uint));

  /* Compute factors */
  N = 0;
  f = 2;
  _A_ = A;

  while (f * f <= _A_) {
    while (_A_ % f == 0) {
      _factors_[N++] = f;
      _A_ /= f;
    }
    f++;
  }
  if (_A_ > 1) _factors_[N++] = _A_;

  *factors = _factors_;
  return 0;
}

// ------------------------------------------------------------------
// INITIALIZE SOURCES
// ------------------------------------------------------------------
int initialize_sources(int Me, int Ntasks, MPI_Comm *Comm, vec2_t mysize,
                       int Nsources, int *Nsources_local, vec2_t **Sources)
{
  srand48(time(NULL) ^ Me);
  int *tasks_with_sources = (int*)malloc(Nsources * sizeof(int));
  
  /* Root task generates which tasks will have sources */
  if (Me == 0) {
    for (int i = 0; i < Nsources; i++)
      tasks_with_sources[i] = (int)lrand48() % Ntasks;
  }
  
  /* Broadcast the assignment */
  MPI_Bcast(tasks_with_sources, Nsources, MPI_INT, 0, *Comm);

  /* Count local sources */
  int nlocal = 0;
  for (int i = 0; i < Nsources; i++)
    nlocal += (tasks_with_sources[i] == Me);
  *Nsources_local = nlocal;
  
  /* Generate local source positions */
  if (nlocal > 0) {
    vec2_t * restrict helper = (vec2_t*)malloc(nlocal * sizeof(vec2_t));      
    for (int s = 0; s < nlocal; s++) {
      helper[s][_x_] = 1 + lrand48() % mysize[_x_];
      helper[s][_y_] = 1 + lrand48() % mysize[_y_];
    }
    *Sources = helper;
  }
  
  free(tasks_with_sources);

  return 0;
}