import numpy as np
import matplotlib.pyplot as plt
import sys

data = np.genfromtxt('DMC.out', skip_header=1)
r = data[:,0]
rho = data[:,1]
eigenvec = data[:,2]
dr = r[1]-r[0]

fig, ax = plt.subplots(figsize=(10,5))

_ = ax.set_ylabel(r'$n(r)$')
_ = ax.set_xlabel(r'$r$')

rho = (rho/r**2.0)/(4.0*np.pi)
condensate = (eigenvec/r)**2.0 / (4.0*np.pi)

_ = ax.plot(r,rho/dr, label=r'$n_{tot}(r)$', marker='o', lw=0.5, alpha=0.75, color='cadetblue')
_ = ax.plot(r,condensate/dr, label=r'$|\phi_0(r)|^2$', marker='o', lw=0.5, alpha=0.75, color='red')

_ = ax.legend()
fig.savefig('Overleaf_obdm.pdf')


data = np.genfromtxt('equilibration_try.out',skip_header=1)

fig, ax = plt.subplots(figsize=(10,5))
steps = data[:,0]
Et = data[:,2]
Nt = data[:,1]
acc_prob = data[:,3]
en_axis = ax.twinx()
_ = ax.plot(steps, Nt, color='firebrick', alpha=0.85, label=r'${Runtime} N_{walk}$')
_ = en_axis.plot(steps,Et, color='cadetblue', alpha=0.75, label=r'${Runtime} \langle E \rangle$')
_ = ax.set_xlabel('Steps')
_ = en_axis.set_ylabel(r'$E_{loc}$', c='cadetblue')
_ = ax.set_ylabel(r'$N_{walk}$', c='firebrick')

#_ = ax.legend() 
#_ = en_axis.legend(loc='center right')

fig.savefig('Overleaf_equilibration.pdf')
