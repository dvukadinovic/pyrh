import matplotlib.pyplot as plt
import numpy as np
import time
import sys

import globin

import pyrh

atmos = globin.Atmosphere("atmos_combined_ss_v2.fits")

start = time.time()


argv = "rhf1d"# -i keyword.input"
argc = len(argv.split(" "))

aux = pyrh.RH(argc, argv)
aux.read_RLK_lines()

for idx in range(1):
	for idy in range(1):
		scale = atmos.data[idx,idy,0]
		temp = atmos.data[idx,idy,1]
		ne = atmos.data[idx,idy,2]
		vz = atmos.data[idx,idy,3]
		vmic = atmos.data[idx,idy,4]
		mag = atmos.data[idx,idy,5]/1e4
		gamma = atmos.data[idx,idy,6]
		chi = atmos.data[idx,idy,7]
		nH = atmos.data[idx,idy,8:]

def run_():
	spec = aux.rhf1d(scale, temp, ne, vz, vmic, mag, gamma, chi, nH, 0)

	return spec

# import timeit
# Nrepeat = 100
# times = timeit.Timer(run_).repeat(1,Nrepeat)
# times = np.array(times)
# print(times)
# print(times/Nrepeat)
# sys.exit()

specRH = run_()

print(time.time() - start)

plt.plot(specRH.I[:-1])
plt.show()

sys.exit()

start = time.time()

run_name = "dummy"
globin.read_input(run_name=run_name)

if globin.mode==0:
	spec = globin.make_synthetic_observations(globin.atm, globin.noise)

print(time.time() - start)

plt.plot(spec.spec[0,0,:,0])
# plt.plot(spec.spec[0,0,:,0] - specRH.I[:-1])

plt.show()