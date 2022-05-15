import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('DMC.out', skip_header=1)
r = data[:,0]
rho = data[:,1]
eigenvec = data[:,2]

fig, ax = plt.subplots()

_ = ax.set_ylabel(r'$n(r)$')
_ = ax.set_xlabel(r'$r$')

rho = (rho/r**2.0)
condensate = (eigenvec/r)**2.0

_ = ax.plot(r,rho, label=r'$n_{tot}(r)$', marker='o', lw=0.5, alpha=0.75)
_ = ax.plot(r,condensate, label=r'$|\phi_0(r)|^2$', marker='o', lw=0.5, alpha=0.75)

plt.legend()
fig.savefig('graph.pdf')
