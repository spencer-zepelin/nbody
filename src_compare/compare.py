import numpy as np
import sys
"""
Tests were performed on the Broadwl partition on Midway with the following commands

128 bodies over 10 timesteps of length 0.1
Parallel was run over 4 ranks with each rank running 2 OMP threads

parallel:
mpiexec -n 4 ./nbody_parallel 128 10 0.1 2

serial:
./nbody_serial 128 10 0.1 1

"""
if __name__ == '__main__':
	assert (len(sys.argv) == 3), "\n---Incorrect Arguments---\nUsage: python3 compare.py <file1> <file2>\n"

	# Test files
	serial_file = sys.argv[1]
	parallel_file = sys.argv[2]

	# Open input files
	fs = open(serial_file, 'rb')
	fp = open(parallel_file, 'rb')

	# Read headers
	first_line = fs.read(8)
	line = np.frombuffer(first_line, dtype=np.int32)
	s_nBodies = int(line[0])
	s_timesteps = int(line[1])

	first_line = fp.read(8)
	line = np.frombuffer(first_line, dtype=np.int32)
	p_nBodies = int(line[0])
	p_timesteps = int(line[1])

	# Ensure headers identical
	assert (s_nBodies == p_nBodies), "Headers differ"
	assert (s_timesteps == p_timesteps), "Headers differ"

	total_entries = s_nBodies * s_timesteps

	acceptable_error = 1.0e-10

	# Loop through each body at every timestep
	for i in range(total_entries):
	    s_line = fs.read(24)
	    p_line = fp.read(24)
	    s_body = np.frombuffer(s_line, dtype=np.float64)
	    p_body = np.frombuffer(s_line, dtype=np.float64)
	    
	    # Assert that serial version of body differs from parallel version by less than acceptable error
	    for j in range(3):
	        assert( abs(s_body[j] - p_body[j]) < acceptable_error), "Data diverges unacceptably"

	# If an assertion error hasn't been thrown, print success message
	print("All divergence within acceptable limits.")
	fs.close()
	fp.close()