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
# 1 -- mass density [cm2/g]
# 2 -- height [km]
atm_scale = int(0)

# wavelength samples for which to compute the spectrum (in nm) in vacuum
wave = np.linspace(630.0, 630.35, num=251) + 0.2

spec = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave)
spec = np.array(spec, dtype=np.float64)

plt.plot(wave, spec[0]*1e8)
plt.xlabel("Wavelength in vacuum [nm]")
plt.ylabel(r"Intensity [10$^{-8}$ W/Hz/srad/m$^2$]")
plt.xlim([wave[0], wave[-1]])
plt.show()