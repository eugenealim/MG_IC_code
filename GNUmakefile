# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_AMRPoisson_execVariableCoefficient

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = Main_PoissonSolver

# the location of the Chombo "lib" directory
#CHOMBO_HOME = /home/cosmos/users/eal40/Chombo2/Chombo-3.2/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTools BoxTools

# the locations of the source code directories
base_dir = .

# input file for 'run' target
INPUT = params.txt

# application-specific targets
src_dirs := Source

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example
