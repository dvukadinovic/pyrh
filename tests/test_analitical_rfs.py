import numpy as np
import matplotlib.pyplot as plt
from time import time
import pyrh

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
    perturbation = 1e-3
    spec_plus = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, loggf_ids=ids, loggf_values=values+perturbation)
    spec_minus = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, loggf_ids=ids, loggf_values=values-perturbation)

    # plt.plot(spec_plus[0], label="plus")
    # plt.plot(spec_minus[0], label="minus")

    rf = (spec_plus[0] - spec_minus[0])/2/perturbation

    print(f"Numerical RFs: {time() - start:.3f}")

    return rf

atmos = np.loadtxt("falc.dat", skiprows=1).T[:,:2]
atmos = spinor2multi(atmos)
scale = atmos[0]

mu = 1.0
cwd = "."
atm_scale = 0 # tau

wave = np.linspace(630.25, 630.5, num=100)

ids = np.array([1], dtype=np.int32)
values = np.array([-0.969], dtype=np.float64)

start = time()
spec = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave, 
					  loggf_ids=ids,
					  loggf_values=values,
                      get_atomic_rfs=True)
print(f"Analytical RFs: {time() - start:.3f}")
spec, rf = spec

rf_num = compute_numerical_rf()

# plt.plot(spec[0], label="orig")

# plt.plot(spec[-1], (-rf[0] - rf_num)/rf_num[0])
plt.plot(spec[-1], rf[0]/spec[0][-1], lw=2, c="C0", label="ana")
ax2 = plt.gca().twinx()
ax2.plot(spec[-1], rf_num/spec[0][-1], lw=0.75, c="k", label="num")
# plt.yscale("log")

plt.legend()

plt.show()