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

start = time.time()

argv = "rhf1d -i keyword.input"
argc = len(argv.split(" "))

scale = atmos.data[0,0,0]
temp = atmos.data[0,0,1]
ne = atmos.data[0,0,2]
vz = atmos.data[0,0,3]
vmic = atmos.data[0,0,4]
mag = atmos.data[0,0,5]
gamma = atmos.data[0,0,6]
chi = atmos.data[0,0,7]
nH = atmos.data[0,0,8:]

j = pyrh.rhf1d(argc, argv, 
			scale, temp, ne, vz, vmic,
			mag, gamma, chi, nH, 0)

print(time.time() - start)

start = time.time()

spec = pyrh.solveray(argc, argv, 
			scale, temp, ne, vz, vmic,
			mag, gamma, chi, nH, 0)

print(time.time() - start)

# plt.plot(spec.I[:-1])

start = time.time()

run_name = "dummy"
globin.read_input(run_name=run_name)

if globin.mode==0:
	spec = globin.make_synthetic_observations(globin.atm, globin.noise)

print(time.time() - start)

# plt.plot(spec.spec[0,0,:,0])

# plt.show()