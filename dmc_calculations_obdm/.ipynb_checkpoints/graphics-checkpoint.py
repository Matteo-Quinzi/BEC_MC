import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('DMC_finer_mesh.out', skip_header=1)
r = data[:,0]
rho = data[:,1]
eigenvec = data[:,2]

fig, ax = plt.subplots(2,1,figsize=(12,12))

_ = ax[0].set_ylabel(r'$n(r)$')
_ = ax[0].set_xlabel(r'$r$')

rho = (rho/r**2.0)
condensate = (eigenvec/r)**2.0

_ = ax[0].plot(r,rho, label=r'$n_{tot}(r)$', marker='o', lw=0.5, alpha=0.75)
_ = ax[0].plot(r,condensate, label=r'$|\phi_0(r)|^2$', marker='o', lw=0.5, alpha=0.75)


data = np.genfromtxt('equilibration_try.out',skip_header=1)
steps = data[:,0]
Et = data[:,2]
Nt = data[:,1]
_ = ax[1].plot(steps, Et)

_= ax[0].legend()
fig.savefig('graph.pdf')
