import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('DMC.out', skip_header=1)
r = data[:,0]
rho = data[:,1]
eigenvec = data[:,2]
dr = r[1]-r[0]

fig, ax = plt.subplots(2,1,figsize=(12,12))

_ = ax[0].set_ylabel(r'$n(r)$')
_ = ax[0].set_xlabel(r'$r$')

rho = (rho/r**2.0)/(4.0*np.pi)
condensate = (eigenvec/r)**2.0 / (4.0*np.pi)

_ = ax[0].plot(r,rho/dr, label=r'$n_{tot}(r)$', marker='o', lw=0.5, alpha=0.75, color='cadetblue')
_ = ax[0].plot(r,condensate/dr, label=r'$|\phi_0(r)|^2$', marker='o', lw=0.5, alpha=0.75, color='red')


data = np.genfromtxt('equilibration_try.out',skip_header=1)
steps = data[:,0]
Et = data[:,2]
Nt = data[:,1]
acc_prob = data[:,3]
walk_axis = ax[1].twinx()
_ = ax[1].plot(steps, Et, color='blue', alpha=0.85, label=r'${Runtime} E_{loc}$')
_ = walk_axis.plot(steps,Nt, color='orange', alpha=0.75, label=r'${Runtime} N_{walk}$')
_ = ax[1].set_xlabel('Steps')
_ = ax[1].set_ylabel(r'$E_{loc}$')
_ = walk_axis.set_ylabel(r'$N_{walk}$')

_ = ax[0].legend()
_ = ax[1].legend() 
_ = walk_axis.legend(loc='center right')
fig.savefig('graph.pdf')
