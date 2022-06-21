# BEC_MC
Project for studying BEC in a gas of N bosons trapped within a harmonic potential.
Hard core interaction between particles is considered.

## VMC 
The VMC procedure allows to study the GS energy using a pair of variational parameters.
It uses a MPI implementation.

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

