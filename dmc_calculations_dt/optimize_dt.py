import numpy as np
import os
import sys

dt = [0.1, 0.01, 0.001, \
      0.0001, 0.00001, 0.000001]

for i in range(len(dt)):
    output_file = 'equilibration_%d.out' % i
    input_file_script = '''
#Number of atoms and walkers and maximum number of walkers
N_at    = 128
N_walk  = 100
N_max   = 100000

#Number of eq.steps / samplings / delta_t between samplings / timestep
eq_it   = 10000
samples = 1000
dt_sam  = 10
dt      = %f
Er      = 1.6883113

#Info about the guiding function
a       = 0.00433
b0      = 0.474
b1      = 0.0014

#Eq.coords input file
coords  = final_coords.txt
eq_out  = %s

''' % (dt[i], output_file)

    #Create temporary input file
    temp_input_file = 'temp_dmc.in'
    os.system('touch %s' % temp_input_file)
    os.system('cat > %s <<EOF%s' % (temp_input_file,input_file_script))

    #Execute dmc run
    os.system('./DMC_omp.x %s %s > log_%d.out' % (temp_input_file, temp_input_file,i))

    #clean input files
    os.system('rm %s' % temp_input_file)
