import numpy as np
import matplotlib.pyplot as plt
from time import time
import pyrh
import sys

def spinor2multi(atm):
	"""
	Casting from SPINOR type atmosphere structure to MULTI type structure.
	"""
	from scipy.constants import k
	new = np.empty(atm.shape, dtype=np.float64)

	# atmosphere scale: height [km] or optical depth @ 500nm
	new[0] = atm[0]
	# temperature [K]
	new[1] = atm[2]
	# electron number dnesity
	new[2] = atm[4]/10/k/atm[2] / 1e6 # [1/cm3]
	# LOS velocity [km/s]
	new[3] = atm[9]/1e5
	# micro-turbulent velocity [km/s]
	new[4] = atm[8]/1e5
	# magnetic field strength [G]
	new[5] = atm[7]
	# magnetic field inclination [rad]
	new[6] = atm[-2]
	# magnetic field azimuth [rad]
	new[7] = atm[-1]
	# total hydrogen number density [1/cm3]
	new[8] = (atm[3] - atm[4])/10/k/atm[2] / 1e6 / 1.26

	return new

def compute_numerical_rf():
	start = time()
	output = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, loggf_ids=ids[:1], loggf_values=values[:1]+perturbation, get_opacities=True)
	spec_plus = output["spectrum"]
	op = output["opacities"]
	output = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, loggf_ids=ids[:1], loggf_values=values[:1]-perturbation, get_opacities=True)
	spec_minus = output["spectrum"]

	# plt.plot(spec_plus[0] - spec_minus[0], label="plus")
	# plt.plot(spec_minus[0], label="minus")
	# plt.show()

	rf = (spec_plus[0] - spec_minus[0])/2/perturbation

	print(f"Numerical RFs: {time() - start:.3f}")

	return rf, op, output["opacities"]

perturbation = 1e-3
atmos = np.loadtxt("falc.dat", skiprows=1).T
atmos = spinor2multi(atmos)
scale = atmos[0]

mu = 1.0
cwd = "."
atm_scale = 0 # tau

wave = np.linspace(630.25, 630.5, num=201)

ids = np.array([0, 1], dtype=np.int32)
values = np.array([-0.71, -0.969], dtype=np.float64)

start = time()

# output = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave,
#                       get_opacities=True
#                     )
# opacities = output["opacities"]
rf_num, op_pos, op_neg = compute_numerical_rf()

output = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, 
					  loggf_ids=ids,
					  loggf_values=values,
                      get_atomic_rfs=True,
                      get_opacities=True
                    )
print(f"Analytical RFs: {time() - start:.3f}")
spec = output["spectrum"]
rf = output["response_functions"]
opacities_anal = output["opacities"]
dopacities_anal = output["d_opacities"]

idw = 150
idp = 0

# print(op_pos[idp,:,idw] - opacities_anal[idp,:,idw])
# print(dopacities_anal[0,idp,:,idw]*perturbation)
# print((op_pos[idp,:,idw] - opacities_anal[idp,:,idw])/dopacities_anal[0,idp,:,idw]/perturbation)

# sys.exit()

# plt.plot(op_pos[0,:,idw] - opacities_anal[0,:,idw], label="diff", c="C1", lw=2)
# plt.plot(dopacities_anal[0,0,:,idw]*perturbation, label="anal", c="C2")
# plt.plot(op_pos[1,:,idw] - opacities_anal[1,:,idw], label="diff", c="C1", lw=2, ls="--")
# plt.plot(dopacities_anal[0,1,:,idw]*perturbation, label="anal", c="C2", ls="--")
# plt.plot(opacities[0,:,idw]*np.log(10)*perturbation, label="theory")
# plt.yscale("log")
# plt.legend()
# plt.show()

# sys.exit()

# print(rf[:,0]/rf_num[:])

# plt.plot(spec[-1], (-rf[0] - rf_num)/rf_num[0])
for idp in range(rf.shape[1]):
	plt.plot(spec[-1], rf[...,idp]/spec[0][-1], lw=2, c="C0", label="ana")
plt.plot(spec[-1], rf_num/spec[0][-1], lw=0.75, c="k", label="num")
# plt.plot(rf[...,0]/rf_num)

# plt.yscale("log")
plt.legend()

plt.show()