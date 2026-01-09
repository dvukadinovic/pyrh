import numpy as np
import matplotlib.pyplot as plt
import sys

nw = 201
nz = 54

B = 1500

# G = np.loadtxt("G_2000.txt").T
# G = G.reshape((4,nw,nz))

# dchi = np.loadtxt("dchi_2000.txt").T
# dchi = dchi.reshape((4,nw,nz))

dI = np.loadtxt(f"dI_{B}G.txt").T
dI = dI.reshape((4,nw,nz))
# dI = np.gradient(dI, axis=2)

Ip = np.loadtxt(f"Ip_{B}G.txt").T
Ip = Ip.reshape((4,nw,nz))
# Ip = np.gradient(Ip, axis=2)

In = np.loadtxt(f"In_{B}G.txt").T
In = In.reshape((4,nw,nz))
# In = np.gradient(In, axis=2)

ids = 3
idz = 2

ana = dI[ids,:,:]
num = (Ip[ids,:,:] - In[ids,:,:])/2/1e-3
continuum = np.max(np.abs(num))
diff = ana - num

# plt.plot(Ip[ids,:,idz], label="Ip")
# plt.plot(In[ids,:,idz], label="In")

# plt.plot(ana, label="analytical")
# plt.plot(num, label="numerical")

# plt.plot(diff/continuum, label="difference")
plt.imshow(diff.T/continuum, aspect='auto', cmap='bwr', vmin=-0.05, vmax=0.05)

plt.legend()

plt.show()