#include "nbody_header.h"

void run_serial_problem(int nBodies, double dt, int nIters, char * fname)
{
	// Open File
	FILE * datafile = fopen(fname, "wb");
	// Check file opened successfully
	assert(datafile != NULL);
	
	// When we open the this binary file for plotting, we will make some assumptions as to
	// size of data types we are writing. As such, we enforce these assumptions here.
	assert(sizeof(int)    == 4 );
	assert(sizeof(double) == 8 );

	// Write Header Info
	fwrite(&nBodies, sizeof(int), 1, datafile);
	fwrite(&nIters,  sizeof(int), 1, datafile);

	// Allocate Bodies
	Body * bodies  = (Body *) calloc( nBodies, sizeof(Body) );
	assert(bodies != NULL);

	// Apply Randomized Initial Conditions to Bodies
	randomizeBodies(bodies, nBodies);

	// Allocate additional space for contiguously stored Cartesian body positions for easier file I/O
	int nPositions = nBodies * 3;
	double * positions = (double *) malloc( nPositions * sizeof(double));
	assert(positions != NULL);

	// Start timer
	double start = get_time();

	// Loop over timesteps
	for (int iter = 0; iter < nIters; iter++)
	{
		printf("iteration: %d\n", iter);

		// Pack up body positions to contiguous buffer
		for( int b = 0, p = 0; b < nBodies; b++ )
		{
			positions[p++] = bodies[b].x;
			positions[p++] = bodies[b].y;
			positions[p++] = bodies[b].z;
		}

		// Output contiguous body positions to file
		fwrite(positions, sizeof(double), nPositions, datafile);

		// Compute new forces & velocities for all particles
		compute_forces(bodies, dt, nBodies);

		// Update positions of all particles
		for (int i = 0 ; i < nBodies; i++)
		{
			bodies[i].x += bodies[i].vx*dt;
			bodies[i].y += bodies[i].vy*dt;
			bodies[i].z += bodies[i].vz*dt;
		}

	}

	// Close data file
	fclose(datafile);

	// Stop timer
	double stop = get_time();

	double runtime = stop-start;
	double time_per_iter = runtime / nIters;
	double interactions = (double) nBodies * (double) nBodies;
	double interactions_per_sec = interactions / time_per_iter;

	printf("SIMULATION COMPLETE\n");
	printf("Runtime [s]:              %.3le\n", runtime);
	printf("Runtime per Timestep [s]: %.3le\n", time_per_iter);
	printf("Interactions per sec:     %.3le\n", interactions_per_sec);

	free(bodies);
	free(positions);
}

void randomizeBodies(Body * bodies, int nBodies)
{
	uint64_t seed = 42;

	// velocity scaling term
	double vm = 1.0e-2;

	for (int i = 0; i < nBodies; i++)
	{
		// Fast forward seed to this particle's location in the global PRNG stream.
		// We forward 6 x particle_id, as each particle requires 7 PRNG samples.
		uint64_t particle_seed = fast_forward_LCG(seed, 7*i);

		// Initialize positions
		bodies[i].x =  2.0 * LCG_random_double(&particle_seed) - 1.0;
		bodies[i].y =  2.0 * LCG_random_double(&particle_seed) - 1.0;
		bodies[i].z = (2.0 * LCG_random_double(&particle_seed) - 1.0) * 0.1;

		// Intialize velocities
		bodies[i].vx =  2.0 * vm * LCG_random_double(&particle_seed) - vm;
		bodies[i].vy =  2.0 * vm * LCG_random_double(&particle_seed) - vm;
		bodies[i].vz = (2.0 * vm * LCG_random_double(&particle_seed) - vm) * 0.1;

		// Give it a spin
		if( bodies[i].x > 0 )
			bodies[i].vy =  5*fabs(bodies[i].vy);
		else
			bodies[i].vy = -5*fabs(bodies[i].vy);

		// Initialize masses so that total mass of system is constant
		// regardless of how many bodies are simulated.
		bodies[i].mass = LCG_random_double(&particle_seed) / nBodies;
	}
}

// Computes the forces between all bodies and updates
// their velocities accordingly
void compute_forces(Body * bodies, double dt, int nBodies)
{
	double G = 6.67259e-3;
	double softening = 1.0e-5;

	// For each particle in the set
	for (int i = 0; i < nBodies; i++)
	{ 
		double Fx = 0.0;
		double Fy = 0.0;
		double Fz = 0.0;

		// Compute force from all other particles in the set
		for (int j = 0; j < nBodies; j++)
		{
			// F_ij = G * [ (m_i * m_j) / distance^3 ] * (location_j - location_i) 

			// First, compute the "location_j - location_i" values for each dimension
			double dx = bodies[j].x - bodies[i].x;
			double dy = bodies[j].y - bodies[i].y;
			double dz = bodies[j].z - bodies[i].z;

			// Then, compute the distance^3 value

			double distance = sqrt(dx*dx + dy*dy + dz*dz + softening);
			double distance_cubed = distance * distance * distance;

			double m_j = bodies[j].mass;
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
