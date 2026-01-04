/*
 * Parallel Stencil Implementation - MPI + OpenMP
 * Adapted to match professor's template structure
 */

#include "stencil_template_parallel.h"

// Global variables
int g_n_omp_threads = 1;
double* g_per_thread_comp_time = NULL;
__thread double thread_local_comp_time = 0.0;

// Function prototypes
int memory_allocate(const int *neighbours, const vec2_t N, buffers_t *buffers_ptr, plane_t *planes_ptr);
int output_energy_stat(int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm);
int initialize(MPI_Comm *Comm, int Me, int Ntasks, int argc, char **argv, vec2_t *S, vec2_t *N,
               int *periodic, int *output_energy_stat_perstep, int *neighbours, int *Niterations,
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
  int Rank, Ntasks;
  int neighbours[4];
  
  int Niterations;
  int periodic;
  vec2_t S, N;
  
  int Nsources;
  int Nsources_local;
  vec2_t *Sources_local;
  double energy_per_source;
  
  plane_t planes[2];
  buffers_t buffers[2];  // buffers[0] = send, buffers[1] = recv
  
  int output_energy_stat_perstep;
  
  /* Initialize MPI */
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
  
  /* Setup OpenMP timing */
  #pragma omp parallel
  {
    #pragma omp master
    { 
      g_n_omp_threads = omp_get_num_threads(); 
    }
  }
  g_per_thread_comp_time = (double*)calloc(g_n_omp_threads, sizeof(double));
  
  /* Initialize */
  int ret = initialize(&myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, 
                       &output_energy_stat_perstep, neighbours, &Niterations,
                       &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
                       &planes[0], &buffers[0]);
  
  if (ret) {
    printf("Task %d terminating with code %d\n", Rank, ret);
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
    double section_start;
    
    /* Inject energy */
    inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
    
    section_start = MPI_Wtime();
    
    /* ====================================== */
    /* HALO EXCHANGE                          */
    /* ====================================== */
    
    double* current_plane = planes[current].data;
    const int sizex = planes[current].size[_x_];
    const int sizey = planes[current].size[_y_];
    const int full_sizex = sizex + 2;
    const int full_sizey = sizey + 2;
    
    // North/South: point directly to data (contiguous)
    buffers[SEND][NORTH] = current_plane + full_sizex;
    buffers[SEND][SOUTH] = current_plane + sizey * full_sizex;
    buffers[RECV][NORTH] = current_plane;
    buffers[RECV][SOUTH] = current_plane + (sizey + 1) * full_sizex;
    
    // East/West: copy to buffers (non-contiguous)
    #pragma omp parallel for
    for (int j = 1; j <= sizey; j++) {
      buffers[SEND][WEST][j-1] = current_plane[j * full_sizex + 1];
      buffers[SEND][EAST][j-1] = current_plane[j * full_sizex + sizex];
    }
    
    /* Non-blocking communication */
    MPI_Request reqs[8];
    int req_idx = 0;
    
    if (neighbours[NORTH] != MPI_PROC_NULL) {
      MPI_Isend(buffers[SEND][NORTH], sizex, MPI_DOUBLE, neighbours[NORTH], 0, myCOMM_WORLD, &reqs[req_idx++]);
      MPI_Irecv(buffers[RECV][NORTH], sizex, MPI_DOUBLE, neighbours[NORTH], 1, myCOMM_WORLD, &reqs[req_idx++]);
    }
    
    if (neighbours[SOUTH] != MPI_PROC_NULL) {
      MPI_Isend(buffers[SEND][SOUTH], sizex, MPI_DOUBLE, neighbours[SOUTH], 1, myCOMM_WORLD, &reqs[req_idx++]);
      MPI_Irecv(buffers[RECV][SOUTH], sizex, MPI_DOUBLE, neighbours[SOUTH], 0, myCOMM_WORLD, &reqs[req_idx++]);
    }
    
    if (neighbours[WEST] != MPI_PROC_NULL) {
      MPI_Isend(buffers[SEND][WEST], sizey, MPI_DOUBLE, neighbours[WEST], 2, myCOMM_WORLD, &reqs[req_idx++]);
      MPI_Irecv(buffers[RECV][WEST], sizey, MPI_DOUBLE, neighbours[WEST], 3, myCOMM_WORLD, &reqs[req_idx++]);
    }
    
    if (neighbours[EAST] != MPI_PROC_NULL) {
      MPI_Isend(buffers[SEND][EAST], sizey, MPI_DOUBLE, neighbours[EAST], 3, myCOMM_WORLD, &reqs[req_idx++]);
      MPI_Irecv(buffers[RECV][EAST], sizey, MPI_DOUBLE, neighbours[EAST], 2, myCOMM_WORLD, &reqs[req_idx++]);
    }
    
    MPI_Waitall(req_idx, reqs, MPI_STATUSES_IGNORE);
    
    total_comm_time += MPI_Wtime() - section_start;
    
    /* Copy East/West halos back */
    #pragma omp parallel for
    for (int j = 1; j <= sizey; j++) {
      current_plane[j * full_sizex] = buffers[RECV][WEST][j-1];
      current_plane[j * full_sizex + sizex + 1] = buffers[RECV][EAST][j-1];
    }
    
    /* Non-periodic boundaries: set to 0 */
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
    /* COMPUTATION                            */
    /* ====================================== */
    
    section_start = MPI_Wtime();
    update_plane(periodic, N, &planes[current], &planes[!current]);
    total_comp_time += MPI_Wtime() - section_start;
    
    if (output_energy_stat_perstep)
      output_energy_stat(iter, &planes[!current], (iter+1) * Nsources * energy_per_source, Rank, &myCOMM_WORLD);
    
    current = !current;
  }
  
  t_end = MPI_Wtime() - t_start;
  
  /* Print timing results */
  double max_comm, max_comp, max_total;
  MPI_Reduce(&total_comm_time, &max_comm, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
  MPI_Reduce(&total_comp_time, &max_comp, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
  MPI_Reduce(&t_end, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
  
  if (Rank == 0) {
    printf("========================================\n");
    printf("Performance Results:\n");
    printf("========================================\n");
    printf("Total time:         %.6f seconds\n", max_total);
    printf("Communication time: %.6f seconds (%.2f%%)\n", max_comm, 100.0 * max_comm / max_total);
    printf("Computation time:   %.6f seconds (%.2f%%)\n", max_comp, 100.0 * max_comp / max_total);
    printf("========================================\n");
  }
  
  output_energy_stat(-1, &planes[current], Niterations * Nsources * energy_per_source, Rank, &myCOMM_WORLD);
  
  /* Cleanup */
  printf("DEBUG: Calling memory_release...\n");
  memory_release(&planes[0], &buffers[0]);
  printf("DEBUG: memory_release done\n");
  
  printf("DEBUG: About to free Sources_local at %p\n", (void*)Sources_local);
  if (Sources_local) free(Sources_local);
  printf("DEBUG: Sources_local freed\n");
  
  printf("DEBUG: About to free g_per_thread_comp_time at %p\n", (void*)g_per_thread_comp_time);
  free(g_per_thread_comp_time);
  printf("DEBUG: g_per_thread_comp_time freed\n");

  printf("DEBUG: Calling MPI_Finalize\n");
  MPI_Finalize();
  return 0;
}


// ------------------------------------------------------------------
// MEMORY ALLOCATE
// ------------------------------------------------------------------
int memory_allocate(const int *neighbours, const vec2_t N, buffers_t *buffers_ptr, plane_t *planes_ptr)
{
  if (!planes_ptr || !buffers_ptr) {
    fprintf(stderr, "Error: Invalid pointers\n");
    return 1;
  }
  
  /* Allocate planes with halo */
  uint frame_size = (planes_ptr[OLD].size[_x_] + 2) * (planes_ptr[OLD].size[_y_] + 2);
  
  planes_ptr[OLD].data = (double*)malloc(frame_size * sizeof(double));
  planes_ptr[NEW].data = (double*)malloc(frame_size * sizeof(double));
  
  if (!planes_ptr[OLD].data || !planes_ptr[NEW].data) {
    perror("malloc failed for planes");
    return 1;
  }
  
  memset(planes_ptr[OLD].data, 0, frame_size * sizeof(double));
  memset(planes_ptr[NEW].data, 0, frame_size * sizeof(double));
  
  /* Allocate East/West buffers */
  int sizey = planes_ptr[OLD].size[_y_];
  
  buffers_ptr[SEND][EAST] = (double*)malloc(sizey * sizeof(double));
  buffers_ptr[SEND][WEST] = (double*)malloc(sizey * sizeof(double));
  buffers_ptr[RECV][EAST] = (double*)malloc(sizey * sizeof(double));
  buffers_ptr[RECV][WEST] = (double*)malloc(sizey * sizeof(double));
  
  if (!buffers_ptr[SEND][EAST] || !buffers_ptr[SEND][WEST] ||
      !buffers_ptr[RECV][EAST] || !buffers_ptr[RECV][WEST]) {
    perror("malloc failed for buffers");
    return 1;
  }
  
  /* North/South point to plane data - set in main loop */
  buffers_ptr[SEND][NORTH] = NULL;
  buffers_ptr[SEND][SOUTH] = NULL;
  buffers_ptr[RECV][NORTH] = NULL;
  buffers_ptr[RECV][SOUTH] = NULL;
  
  return 0;
}

// ------------------------------------------------------------------
// MEMORY RELEASE
// ------------------------------------------------------------------
int memory_release(plane_t *planes, buffers_t *buffers)
{
  if (planes) {
    if (planes[OLD].data) free(planes[OLD].data);
    if (planes[NEW].data) free(planes[NEW].data);
  }
  
  /* Free ONLY East/West buffers (North/South are pointers to plane data) */
  if (buffers) {
    printf("DEBUG: About to free EAST buffers...\n");
    if (buffers[SEND][EAST]) {
      printf("DEBUG: Freeing SEND EAST at %p\n", buffers[SEND][EAST]);
      free(buffers[SEND][EAST]);
    }
    if (buffers[SEND][WEST]) {
      printf("DEBUG: Freeing SEND WEST at %p\n", buffers[SEND][WEST]);
      free(buffers[SEND][WEST]);
    }
    printf("DEBUG: About to free RECV buffers...\n");
    if (buffers[RECV][EAST]) {
      printf("DEBUG: Freeing RECV EAST at %p\n", buffers[RECV][EAST]);
      free(buffers[RECV][EAST]);
    }
    if (buffers[RECV][WEST]) {
      printf("DEBUG: Freeing RECV WEST at %p\n", buffers[RECV][WEST]);
      free(buffers[RECV][WEST]);
    }
    printf("DEBUG: Buffer cleanup done\n");
  }
  
  return 0;
}



// ------------------------------------------------------------------
// OUTPUT ENERGY STAT
// ------------------------------------------------------------------
int output_energy_stat(int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm)
{
  double local_energy = 0, total_energy = 0;
  get_total_energy(plane, &local_energy);
  
  MPI_Reduce(&local_energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm);
  
  if (Me == 0) {
    if (step >= 0)
      printf("[ step %4d ] ", step);
    printf("injected: %g, system: %g (avg: %g per point)\n",
           budget, total_energy, total_energy / (plane->size[_x_] * plane->size[_y_]));
  }
  
  return 0;
}

// ------------------------------------------------------------------
// INITIALIZE
// ------------------------------------------------------------------
int initialize(MPI_Comm *Comm, int Me, int Ntasks, int argc, char **argv,
               vec2_t *S, vec2_t *N, int *periodic, int *output_energy_stat_perstep,
               int *neighbours, int *Niterations, int *Nsources, int *Nsources_local,
               vec2_t **Sources_local, double *energy_per_source,
               plane_t *planes, buffers_t *buffers)
{
  int verbose = 0;
  vec2_t Grid;
  
  /* Default values */
  (*S)[_x_] = 1000;
  (*S)[_y_] = 1000;
  *Nsources = 1;
  *energy_per_source = 1.0;
  *Niterations = 100;
  *periodic = 0;
  *output_energy_stat_perstep = 0;
  
  /* Parse arguments */
  int opt;
  while ((opt = getopt(argc, argv, "x:y:e:E:n:p:o:v:h")) != -1) {
    switch (opt) {
      case 'x': (*S)[_x_] = atoi(optarg); break;
      case 'y': (*S)[_y_] = atoi(optarg); break;
      case 'e': *Nsources = atoi(optarg); break;
      case 'E': *energy_per_source = atof(optarg); break;
      case 'n': *Niterations = atoi(optarg); break;
      case 'p': *periodic = atoi(optarg); break;
      case 'o': *output_energy_stat_perstep = atoi(optarg); break;
      case 'v': verbose = atoi(optarg); break;
      case 'h':
        if (Me == 0) {
          printf("Usage: %s [options]\n"
                 "  -x <size>  x size [1000]\n"
                 "  -y <size>  y size [1000]\n"
                 "  -e <num>   energy sources [1]\n"
                 "  -E <val>   energy per source [1.0]\n"
                 "  -n <iter>  iterations [100]\n"
                 "  -p <0|1>   periodic boundaries [0]\n"
                 "  -o <0|1>   output every step [0]\n"
                 "  -v <level> verbosity [0]\n", argv[0]);
        }
        MPI_Finalize();
        exit(0);
    }
  }
  
  /* Validate */
  if ((*S)[_x_] <= 0 || (*S)[_y_] <= 0 || *Nsources <= 0 || *Niterations <= 0) {
    if (Me == 0) printf("Error: Invalid parameters\n");
    return 1;
  }
  
  /* Domain decomposition - make grid as square as possible */
  Grid[_x_] = 1;
  Grid[_y_] = Ntasks;
  int sqrt_n = (int)sqrt(Ntasks);
  for (int i = sqrt_n; i > 0; i--) {
    if (Ntasks % i == 0) {
      Grid[_x_] = i;
      Grid[_y_] = Ntasks / i;
      break;
    }
  }
  
  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];
  
  /* Task coordinates */
  int X = Me % Grid[_x_];
  int Y = Me / Grid[_x_];
  
  /* Set neighbours */
  neighbours[NORTH] = (Y > 0) ? Me - Grid[_x_] : 
                      (*periodic ? Me + (Grid[_y_] - 1) * Grid[_x_] : MPI_PROC_NULL);
  neighbours[SOUTH] = (Y < Grid[_y_] - 1) ? Me + Grid[_x_] : 
                      (*periodic ? Me - (Grid[_y_] - 1) * Grid[_x_] : MPI_PROC_NULL);
  neighbours[WEST] = (X > 0) ? Me - 1 : 
                     (*periodic ? Me + (Grid[_x_] - 1) : MPI_PROC_NULL);
  neighbours[EAST] = (X < Grid[_x_] - 1) ? Me + 1 : 
                     (*periodic ? Me - (Grid[_x_] - 1) : MPI_PROC_NULL);
  
  /* Compute local size */
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
  
  if (verbose && Me == 0) {
    printf("Grid: %d x %d tasks, Global size: %d x %d\n", 
           Grid[_x_], Grid[_y_], (*S)[_x_], (*S)[_y_]);
  }
  
  /* Allocate memory */
  int ret = memory_allocate(neighbours, *N, buffers, planes);
  if (ret) return ret;
  
  /* Initialize sources */
  ret = initialize_sources(Me, Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local);
  return ret;
}

// ------------------------------------------------------------------
// SIMPLE FACTORIZATION
// ------------------------------------------------------------------
uint simple_factorization(uint A, int *Nfactors, uint **factors)
{
  int N = 0;
  uint _A_ = A;
  
  for (int f = 2; f * f <= _A_; f++) {
    while (_A_ % f == 0) {
      N++;
      _A_ /= f;
    }
  }
  if (_A_ > 1) N++;
  
  *Nfactors = N;
  uint *_factors_ = (uint*)malloc(N * sizeof(uint));
  
  N = 0;
  _A_ = A;
  for (int f = 2; f * f <= _A_; f++) {
    while (_A_ % f == 0) {
      _factors_[N++] = f;
      _A_ /= f;
    }
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
  int *tasks = (int*)malloc(Nsources * sizeof(int));
  
  if (Me == 0) {
    for (int i = 0; i < Nsources; i++)
      tasks[i] = lrand48() % Ntasks;
  }
  
  MPI_Bcast(tasks, Nsources, MPI_INT, 0, *Comm);
  
  int nlocal = 0;
  for (int i = 0; i < Nsources; i++)
    nlocal += (tasks[i] == Me);
  
  *Nsources_local = nlocal;
  
  if (nlocal > 0) {
    vec2_t *helper = (vec2_t*)malloc(nlocal * sizeof(vec2_t));
    for (int s = 0; s < nlocal; s++) {
      helper[s][_x_] = 1 + lrand48() % mysize[_x_];
      helper[s][_y_] = 1 + lrand48() % mysize[_y_];
    }
    *Sources = helper;
  }
  
  free(tasks);
  return 0;
}
