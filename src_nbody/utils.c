#include "nbody_header.h"
double get_time(void)
{
	#ifdef MPI
	return MPI_Wtime();
	#endif

	#ifdef OPENMP
	return omp_get_wtime();
	#endif

	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	long ms = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
	double time = (double) ms / 1000.0;

	return time;
}

void print_inputs(long nBodies, double dt, int nIters, int nthreads )
{
	int mype = 0;
	int nprocs = 1;

	#ifdef MPI
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif

	if( mype == 0 )
	{
		printf("INPUT PARAMETERS:\n");
		printf("N Bodies =                     %ld\n", nBodies);
		printf("Timestep dt =                  %.3le\n", dt);
		printf("Number of Timesteps =          %d\n", nIters);
		printf("Number of Threads per Rank =   %d\n", nthreads);
		#ifdef MPI
		printf("Number of MPI Ranks =          %d\n", nprocs);
		#endif
		printf("BEGINNING N-BODY SIMLUATION\n");
	}
}

// A 63-bit LCG
// Returns a double precision value from a uniform distribution
// between 0.0 and 1.0 using a caller-owned state variable.
double LCG_random_double(uint64_t * seed)
{
	const uint64_t m = 9223372036854775808ULL; // 2^63
	const uint64_t a = 2806196910506780709ULL;
	const uint64_t c = 1ULL;
	*seed = (a * (*seed) + c) % m;
	return (double) (*seed) / (double) m;
}

// "Fast Forwards" an LCG PRNG stream
// seed: starting seed
// n: number of iterations (samples) to forward
// Returns: forwarded seed value
uint64_t fast_forward_LCG(uint64_t seed, uint64_t n)
{
	const uint64_t m = 9223372036854775808ULL; // 2^63
	uint64_t a = 2806196910506780709ULL;
	uint64_t c = 1ULL;

	n = n % m;

	uint64_t a_new = 1;
	uint64_t c_new = 0;

	while(n > 0)
	{
		if(n & 1)
		{
			a_new *= a;
			c_new = c_new * a + c;
		}
		c *= (a + 1);
		a *= a;

		n >>= 1;
	}
	return (a_new * seed + c_new) % m;
}
