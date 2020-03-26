#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
OPTIMIZE    = yes
MPI         = yes
OPENMP      = yes
DEBUG       = no
PROFILE     = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = nbody

source = \
main.c \
utils.c \
serial.c \
parallel.c 

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# MPI Compiler
ifeq ($(MPI),yes)
  CC = mpicc
  CFLAGS += -DMPI
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg -fno-omit-frame-pointer
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

# OpenMP
ifeq ($(OPENMP),yes)
  CFLAGS += -fopenmp -DOPENMP
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) nbody_header.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c Makefile nbody_header.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) *\.lst nbody.dSYM

edit:
	vim -p $(source) nbody_header.h

run: $(program)
	./$(program)

mpi: $(program)
	mpirun -n 2 ./$(program)
