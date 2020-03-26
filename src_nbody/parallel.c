#include "nbody_header.h"

#ifdef MPI
void run_parallel_problem(int nBodies, double dt, int nIters, char * fname)
{
	// MPI initialization
	int nprocs; // number of processes (aka "procs") used in this invocation
    int mype; // my process id (from 0 .. nprocs-1)
    int stat; // used as error code for MPI calls
    MPI_Status status; // MPI specific error variable
    stat = MPI_Comm_size(MPI_COMM_WORLD, &nprocs); // get number of procs
    assert(stat == MPI_SUCCESS); // Check to make sure call worked
    stat = MPI_Comm_rank(MPI_COMM_WORLD, &mype); // get my integer proc id
    assert(stat == MPI_SUCCESS); // Check to make sure call worked

    // Cartesian ring
    MPI_Comm ring_comm;
    int true=1, left, right;
    MPI_Cart_create( MPI_COMM_WORLD, 1, &nprocs, &true, 1, &ring_comm );
  	MPI_Cart_shift( ring_comm, 0, 1, &left, &right );

	// Open File
	MPI_File datafile;
	MPI_File_open( MPI_COMM_WORLD, fname, 
	    	MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &datafile);
	// Check file opened successfully
	assert(datafile != NULL);
	
	// When we open the this binary file for plotting, we will make some assumptions as to
	// size of data types we are writing. As such, we enforce these assumptions here.
	assert(sizeof(int)    == 4 );
	assert(sizeof(double) == 8 );

	// Write Header Info
	// Only rank zero!
	if (mype == 0){
		MPI_File_write(datafile, &nBodies, 1, MPI_INT, &status);
		MPI_File_write(datafile, &nIters, 1, MPI_INT, &status);
	}

	// Assume even distribution of bodies over ranks
	assert(nBodies % nprocs == 0);
	int nBodies_per_rank = nBodies / nprocs;

	// Allocate Bodies
	Body * bodies  = (Body *) calloc( nBodies_per_rank, sizeof(Body) );
	assert(bodies != NULL);

	// Apply Randomized Initial Conditions to Bodies
	parallel_randomizeBodies(bodies, nBodies, nBodies_per_rank, mype, nprocs);


	// NOTE:
	// We only need to write position, but we need to send position and mass

	// Allocate additional space for contiguously stored Cartesian body positions for easier file I/O
	int nPositions_per_rank = nBodies_per_rank * 3;
	double * positions = (double *) malloc( nPositions_per_rank * sizeof(double));
	assert(positions != NULL);

	// Allocate space for send data
	int nPositionmass_per_rank = nBodies_per_rank * 4;
	double * send_buf = (double *) malloc( nPositionmass_per_rank * sizeof(double));
	double * recv_buf = (double *) malloc( nPositionmass_per_rank * sizeof(double));
	assert(send_buf != NULL && recv_buf != NULL);

	// Start timer
	double start = get_time();

/* MAIN LOOP */ 
	for (int iter = 0; iter < nIters; iter++)
	{
		if (mype == 0){
			printf("iteration: %d\n", iter);
		}
		
		// Pack up body positions to write buffer and send buffer
#ifdef OPENMP
#pragma omp parallel for shared(positions, send_buf, bodies) schedule(static)
#endif
		for( int b = 0; b < nBodies_per_rank; b++ )
		{
			int p = b * 3;
			int r = b * 4;
			positions[p++] = bodies[b].x;
			positions[p++] = bodies[b].y;
			positions[p] = bodies[b].z;
			send_buf[r++] = bodies[b].x;
			send_buf[r++] = bodies[b].y;
			send_buf[r++] = bodies[b].z;
			send_buf[r] = bodies[b].mass;
		}

		// Collectively write body positions to file
		distributed_write_timestep(positions, nBodies, nBodies_per_rank, iter, mype, &datafile, status);


		// Perform force/velocity calc of own bodies
		compute_forces_multi_set(bodies, send_buf, dt, nBodies_per_rank, 1);

		// Pipeline
		for (int push = 0; push < nprocs-1; push++){
			// Send left; recv from right
			MPI_Sendrecv(send_buf, nPositionmass_per_rank, MPI_DOUBLE, left, 99, 
						recv_buf, nPositionmass_per_rank, MPI_DOUBLE, right, MPI_ANY_TAG, ring_comm, &status);
			int count;
			MPI_Get_count(&status, MPI_DOUBLE, &count);
			// Shift memory
			memcpy(send_buf, recv_buf, nPositionmass_per_rank * sizeof(double));
			// Perform force/velocity calc on new data
			compute_forces_multi_set(bodies, send_buf, dt, nBodies_per_rank, 0);
		}

		// Update positions of all particles
#ifdef OPENMP
#pragma omp parallel for shared(bodies) schedule(static)
#endif
		for (int i = 0 ; i < nBodies_per_rank; i++)
		{
			bodies[i].x += bodies[i].vx*dt;
			bodies[i].y += bodies[i].vy*dt;
			bodies[i].z += bodies[i].vz*dt;
		}

	}

	// Close data file
	MPI_File_close(&datafile);

	// Stop timer
	double stop = get_time();

	double runtime = stop-start;
	double time_per_iter = runtime / nIters;
	double interactions = (double) nBodies * (double) nBodies;
	double interactions_per_sec = interactions / time_per_iter;

	if (mype == 0){
		printf("SIMULATION COMPLETE\n");
		printf("Runtime [s]:              %.3le\n", runtime);
		printf("Runtime per Timestep [s]: %.3le\n", time_per_iter);
		printf("Interactions per sec:     %.3le\n", interactions_per_sec);		
	}

	free(bodies);
	free(positions);
	free(send_buf);
	free(recv_buf);
}

void compute_forces_multi_set(Body * bodies, double * remote, double dt, int nBodies_per_rank, int self)
{
	double G = 6.67259e-3;
	double softening = 1.0e-5;

	// For each particle in the local
#ifdef OPENMP
#pragma omp parallel for shared(remote, bodies) schedule(static)
#endif
	for (int i = 0; i < nBodies_per_rank; i++)
	{ 
		double Fx = 0.0;
		double Fy = 0.0;
		double Fz = 0.0;

		// Compute force from all other particles in the remote
		for (int j = 0; j < nBodies_per_rank; j++)
		{
			// Unless computing force on self
			if (self && i == j){
				continue;
			}

			double dx = remote[4 * j] - bodies[i].x;
			double dy = remote[(4* j) + 1] - bodies[i].y;
			double dz = remote[(4 * j) + 2] - bodies[i].z;

			// distance = sqrt( dx^2 + dx^2 + dz^2 )
			double distance = sqrt(dx*dx + dy*dy + dz*dz + softening);
			double distance_cubed = distance * distance * distance;

			double m_j = remote[(4 * j) + 3];
			double mGd = G * m_j / distance_cubed;
			Fx += mGd * dx;
			Fy += mGd * dy;
			Fz += mGd * dz;
		}

		bodies[i].vx += dt*Fx;
		bodies[i].vy += dt*Fy;
		bodies[i].vz += dt*Fz;
	}
}



void parallel_randomizeBodies(Body * bodies, int nBodies, int nBodies_per_rank, int mype, int nprocs)
{
	uint64_t seed = 42;

	// velocity scaling term
	double vm = 1.0e-2;
#ifdef OPENMP
#pragma omp parallel for shared(bodies) schedule(static)
#endif
	for (int i = 0; i < nBodies_per_rank; i++)
	{
		int global_particle_id = (mype * nBodies_per_rank) + i;
		// Fast forward seed to this particle's location in the global PRNG stream.
		// We forward 4 x particle_id, as each particle requires 4 PRNG samples.
		uint64_t particle_seed = fast_forward_LCG(seed, global_particle_id * 4);

		// Initialize positions
		double d1 = LCG_random_double(&particle_seed);
		double d2 = LCG_random_double(&particle_seed);
		double d3 = LCG_random_double(&particle_seed);

		bodies[i].x =  2.0 * d1 - 1.0;
		bodies[i].y =  2.0 * d2 - 1.0;
		bodies[i].z = (2.0 * d3 - 1.0) * 0.1;

		// Intialize velocities
		bodies[i].vx =  3 * (2.0 * vm * d1 - vm);
		bodies[i].vy =  3 * (2.0 * vm * d2 - vm);
		bodies[i].vz = (2.0 * vm * d3 - vm) * 0.1;


		// Initialize masses so that total mass of system is constant
		// regardless of how many bodies are simulated.
		bodies[i].mass = LCG_random_double(&particle_seed) / nBodies;
	}
}
// Writes all particle locations for a single timestep
void distributed_write_timestep(double * positions, int nBodies, int nBodies_per_rank, int timestep, int mype, MPI_File * fh, MPI_Status status)
{
	int header = 1;
	int bytes_per_step = nBodies * 3; 
	int bytes_per_rank = nBodies_per_rank * 3; 
	// Offset is number of doubles
	MPI_Offset offset = sizeof(double) * (header + (bytes_per_step * timestep) + (mype * bytes_per_rank));
	// Collective write
	MPI_File_write_at_all(*fh, offset, positions, nBodies_per_rank * 3, MPI_DOUBLE, &status);

}
#endif
