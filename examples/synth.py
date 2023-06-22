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

atmos = np.loadtxt("falc.dat", skiprows=1)
atmos = np.array(atmos.T, dtype=np.float64)
atmos= spinor2multi(atmos)

# path to the *.input files
cwd = "."

# mu angle for which to compute the spectrum
mu = 1.0

# type of atmosphere stratification:
# 0 -- optical depth @ 500nm
# 1 -- mass density
# 2 -- height [km]
atm_scale = int(0)

# wavelength samples for which to compute the spectrum (in nm) in vacuum
wave = np.linspace(630.0, 630.35, num=251) + 0.2

spec = pyrh.compute1d(cwd, mu, atm_scale, atmos,
					  wave,
					  int(0), np.empty(1, dtype=np.float64), np.empty((3,1), dtype=np.float64),
					  np.array([], dtype=np.int32), np.array([], dtype=np.float64),
					  np.array([], dtype=np.int32), np.array([], dtype=np.float64))
spec = np.array(spec, dtype=np.float64)

plt.plot(wave, spec[0])
plt.show()