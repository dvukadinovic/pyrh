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
	pyrh.compute1d(cwd, mu, atm_scale, spec_plus, atmos, wave, loggf_ids=IDs[:1], loggf_values=values[:1]+perturbation)
	pyrh.compute1d(cwd, mu, atm_scale, spec_minus, atmos, wave, loggf_ids=IDs[:1], loggf_values=values[:1]-perturbation)

	# plt.plot(spec_plus[:,0] - spec_minus[:,0], label="diff")
	# plt.plot(spec_plus[:,0], label="plus")
	# plt.plot(spec_minus[:,0], label="minus")
	# plt.show()

	rf = np.empty((4, spec_plus[:,0].shape[0]), dtype=np.float64)
	for _ids in range(4):
		rf[_ids] = (spec_plus[:,_ids] - spec_minus[:,_ids])/2/perturbation

	print(f"Numerical RFs: {time() - start:.3f}")

	return rf

perturbation = 1e-3
atmos = np.loadtxt("falc.dat", skiprows=1).T
atmos = spinor2multi(atmos)

atmos[6] = np.deg2rad(45.0)  # inclination
atmos[7] = np.deg2rad(10.0)  # azimuth
atmos[4] = 0.0  # microturbulence

mu = 1.0
cwd = "."
atm_scale = 0 # tau

wave = np.linspace(630.255, 630.5, num=201)

IDs = np.array([0], dtype=np.int32)
values = np.array([-0.71], dtype=np.float64)

spec = np.zeros((len(wave), 4), dtype=np.float64)
spec_plus = np.zeros((len(wave), 4), dtype=np.float64)
spec_minus = np.zeros((len(wave), 4), dtype=np.float64)

rfs_analytical = np.zeros((len(wave), 4, len(values)), dtype=np.float64)
analytical
ids = 0

# for idb, B in enumerate([0,500,2000,8000]):
for idb, B in enumerate([1500]):
	atmos[5] = B

	scale = atmos[0]

	rf_num = compute_numerical_rf()
	print("----")

	start = time()
	out = pyrh.compute1d(cwd, mu, atm_scale, spec, atmos, wave, 
						loggf_ids=IDs,
						loggf_values=values,
						get_atomic_rfs=True,
						rfs=rfs_analytical
						)
	print(f"Analytical RFs: {time() - start:.3f}")

	print("==========================")

	plt.plot(rfs_analytical[:,ids,0]/spec[0][0], c=f"C{idb}", label="Analytical RF")
	plt.plot(rf_num[ids]/spec[0][0], linestyle="--", c=f"C{idb}", label="Numerical RF")

	# plt.plot((rfs_analytical[0,:,ids]/rf_num[ids] - 1)*100, label=f"rel. diff B={B}")
	# maxRf = np.max(np.abs(rf_num[0]))
	# plt.plot((rfs_analytical[0,:,ids] - rf_num[ids])/maxRf, label=B)

	# plt.yscale("symlog")

	print(rf_num[:,49], rfs_analytical[49,0,:])

plt.legend()

plt.show()