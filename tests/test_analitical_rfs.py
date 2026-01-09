import numpy as np
import matplotlib.pyplot as plt
from time import time
import pyrh
import sys

# stokes = np.load("stokes_spec.npy")
# I = np.load("nonpol_spec.npy")

# print(stokes[59]/stokes[0], I[59]/I[0])

# plt.plot(stokes, label="pol")
# plt.plot(I, label="non pol")
# plt.legend()
# plt.show()

# sys.exit()

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
	# spec_plus = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, loggf_ids=ids[:1], loggf_values=values[:1]+perturbation)
	spec_minus = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, loggf_ids=ids[:1], loggf_values=values[:1]-perturbation)

	# plt.plot(spec_plus[0] - spec_minus[0], label="plus")
	# plt.plot(spec_minus[0], label="minus")
	# plt.show()

	rf = np.empty((4, spec_plus[0].shape[0]), dtype=np.float64)
	for _ids in range(4):
		rf[_ids] = (spec_plus[_ids] - spec_minus[_ids])/2/perturbation

	print(f"Numerical RFs: {time() - start:.3f}")

	return rf

perturbation = 1e-3
atmos = np.loadtxt("falc.dat", skiprows=1).T
atmos = spinor2multi(atmos)

atmos[6] = np.deg2rad(45.0)  # inclination
atmos[7] = np.deg2rad(10.0)   # az
# atmos[3] = (np.arange(atmos.shape[1]) - atmos.shape[1]/2)/15  # LOS vel
atmos[4] = 0.0  # microturbulence

# for idb, B in enumerate([0,500,2000,8000]):
for idb, B in enumerate([1500]):
	atmos[5] = B

	scale = atmos[0]

	mu = 1.0
	cwd = "."
	atm_scale = 0 # tau

	wave = np.linspace(630.255, 630.5, num=201)

	ids = np.array([0], dtype=np.int32)
	values = np.array([-0.71], dtype=np.float64)

	rf_num = compute_numerical_rf()
	print("----")

	start = time()
	spec, rf = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, 
						loggf_ids=ids,
						loggf_values=values,
						get_atomic_rfs=True
						)
	print(f"Analytical RFs: {time() - start:.3f}")

	print("==========================")
	continue

	ids = 0

	plt.plot(rf[0,:,ids]/spec[0][0], c=f"C{idb}", label="Analytical RF")
	plt.plot(rf_num[ids]/spec[0][0], linestyle="--", c=f"C{idb}", label="Numerical RF")

	# plt.plot((rf[0,:,ids]/rf_num[ids] - 1)*100, label=f"rel. diff B={B}")
	# maxRf = np.max(np.abs(rf_num[0]))
	# plt.plot((rf[0,:,ids] - rf_num[ids])/maxRf, label=B)

	# plt.yscale("symlog")

	print(rf_num[:,49], rf[0,49,:])

plt.legend()

plt.show()