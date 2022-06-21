# BEC_MC
Project for studying BEC in a gas of N bosons trapped within a harmonic potential.\
Hard core interaction between particles is considered.

## How to use this repository
- `code` contains the source code + some makefile to compile the source using ifort
- `exec_files` contains the executables already compiled (using `ifort`) + a subdirectory with examples of input/outpuf file and `.log` files showing what the programs give as diagnostic informations
- Other directories contain intermediate results (nothing too useful in general, but you can have a look if curious)
- `report.pdf` contains a detailed description of what is done
- To have the executable files work correctly please respect the structure used in the example input files and please have all input file needed in the same directory of the executable file.

Here follows a description of the input/output files (if needed for their interpretation).

## VMC 
The VMC procedure allows to study the GS energy using a pair of variational parameters.
It uses a MPI implementation (Fized number of walkers).
Example of usage:
`mpirun -np 4 ./VMC.x VMC.in VMC.out`

### VMC input file
N = number of bosons
N = number of Monte Carlo walkers
M = number of steps per each cycle
N = number of cycles
a = ratio of the scattering length to trap length
b = first variational parameter
b = second variational parameter
gas_radius = maximum initialization radius for the gas 
from_file = boolean telling if initial configuration should be taken from a file
coords_file = if from_file == 1 it gives the name of the file with the initial coordinates

### VMC outputs
- Starting configuration of the gas (same for all walkers)
- Final configuration of one walker 
- Runtime local energy of the last cycle for each mpi process
- File with the average value of the energy, average error and average acceptance ratio on each cycle 

## DMC
The DMC procedure explores the configuration space using a guiding function.
It uses an OpenMP implementation (variable number of walkers).
Example of usage:
`./DMC.x DMC.in DMC.out`

### DMC input file
N_at = number of bosons
N_walk = initial population of walkers
N_max = maximum population allowed of walkers
eq_it = number of iterations for the equilibration phase
samples = number of samples of the many-body wavefunction in the sampling phase
dt_sam = numper of steps between two successive samples
dt = timestep
Er = initial reference value for the energy scale
a = ratio of scattering length to trap length
b0 = first variational parameter
b1 = second variational parameter
coords = file storing a initial configuration of coordinates (taken from a VMC procedure)
eq_out = name of the output file with the info about the equilibration phase
Nl = number of points in the discretization of the radial distance
r_min = minimum distance from the center of the trap
r_max = maximum distance from the center of the trap


### DMC output files
- An output file containing the values of the radial distance and particle density and condensate orbital evaluated at those values
- A file containing the values of the population of walkers, energy and acceptance probability for each step of the equilibration procedure
