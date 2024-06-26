import numpy as np
import matplotlib.pyplot as plt

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

atmos = np.loadtxt("falc.dat", skiprows=1).T
atmos = spinor2multi(atmos)

scale = atmos[0]

cwd = "."
atm_scale = 0 # tau

old = atmos[2].copy()

wave = np.linspace(630.25, 630.5, num=1001)

spec = pyrh.compute1d(cwd, 1.0, atm_scale, atmos, wave)

plt.plot(spec[-1], spec[0])
plt.show()