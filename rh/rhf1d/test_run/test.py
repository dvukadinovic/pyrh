import matplotlib.pyplot as plt
import numpy as np
import time
import sys

import globin

import pyrh

# argv = "rhf1d -i keyword.input"
# argc = len(argv.split(" "))

# inData = pyrh.read_input(argc, argv)
# print()
# # print(inData["magneto_optical"])
# print(inData)

# sys.exit()

atmos = globin.Atmosphere("atm_0_0.atmos")
atmos = globin.Atmosphere("atmos_combined_ss_v2.fits")

start = time.time()

argv = "rhf1d -i keyword.input"
argc = len(argv.split(" "))

for idx in range(atmos.nx):
	for idy in range(atmos.ny):
		scale = atmos.data[idx,idy,0]
		temp = atmos.data[idx,idy,1]
		ne = atmos.data[idx,idy,2]
		vz = atmos.data[idx,idy,3]
		vmic = atmos.data[idx,idy,4]
		mag = atmos.data[idx,idy,5]
		gamma = atmos.data[idx,idy,6]
		chi = atmos.data[idx,idy,7]
		nH = atmos.data[idx,idy,8:]

		j = pyrh.rhf1d(argc, argv, 
					scale, temp, ne, vz, vmic,
					mag, gamma, chi, nH, 0)

		spec = pyrh.solveray(argc, argv, 
					scale, temp, ne, vz, vmic,
					mag, gamma, chi, nH, 0)
# plt.plot(spec.I[:-1])
# plt.show()

print(time.time() - start)


start = time.time()

run_name = "dummy"
globin.read_input(run_name=run_name)

if globin.mode==0:
	spec = globin.make_synthetic_observations(globin.atm, globin.noise)

print(time.time() - start)

# plt.plot(spec.spec[0,0,:,0])

# plt.show()