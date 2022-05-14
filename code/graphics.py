import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('DMC.out', skip_header=1)
r = data[:,0]
rho = data[:,1]
eigenvec = data[:,2]

fig, ax = plt.subplots()
_ = ax.plot(r,rho/max(rho)/r**2.0)
_ = ax.plot(r,eigenvec/max(eigenvec)/r)
fig.savefig('graph.pdf')
