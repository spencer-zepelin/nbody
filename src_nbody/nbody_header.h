#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#ifdef MPI
#include <mpi.h>
#include <memory.h>
#endif
#ifdef OPENMP
#include <omp.h>
#endif

typedef struct {
	// Location r_i = (x,y,z)
	double x;
	double y;
	double z;
	// Velocity v_i = (vx, vy, vz)
	double vx;
	double vy;
	double vz;
	// Mass
	double mass;
} Body;

// serial.c
void run_serial_problem(int nBodies, double dt, int nIters, char * fname);
void randomizeBodies(Body * bodies, int nBodies);
void randomizeBodies_fun(Body * bodies, int nBodies);
void compute_forces(Body * bodies, double dt, int nBodies);

// utils.c
double get_time(void);
void print_inputs(long nBodies, double dt, int nIters, int nthreads );
double LCG_random_double(uint64_t * seed);
uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);

// parallel.c
#ifdef MPI
void compute_forces_multi_set(Body * bodies, double * remote, double dt, int nBodies_per_rank, int self);
void parallel_randomizeBodies(Body * bodies, int nBodies, int nBodies_per_rank, int mype, int nprocs);
void distributed_write_timestep(double * positions, int nBodies, int nBodies_per_rank, int timestep, int mype, MPI_File * fh, MPI_Status status);
void run_parallel_problem(int nBodies, double dt, int nIters, char * fname);
#endif
