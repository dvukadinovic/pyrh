import matplotlib.pyplot as plt
import numpy as np
import pyrh
import sys

argv = "rhf1d -i keyword.input"
argc = len(argv.split(" "))

inData = pyrh.read_input(argc, argv)
print()
# print(inData["magneto_optical"])
print(inData)

sys.exit()

import time

import globin
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

spec = pyrh.py_rhf1d(argc, argv, 
			scale, temp, ne, vz, vmic,
			mag, gamma, chi, nH, 0)

print(time.time() - start)
# plt.plot(spec.lam, spec.I[:,-1])
# plt.show()