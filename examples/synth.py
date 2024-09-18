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
	new[8] = (atm[3] - atm[4])/10/k/atm[2] / 1e6 / 1.24

	return new

def air_to_vacuum(wavelength, air2vacuum_limit=199.9352):
    sigma_sq = (1.0e7/wavelength)**2
    fact = 1.0000834213 + 2.406030e6/(1.3e10 - sigma_sq) + 1.5997e4/(3.89e9 - sigma_sq)

    ind = np.argmin(abs(wavelength-air2vacuum_limit))
    fact[wavelength<air2vacuum_limit] = 1

    aux = wavelength*fact

    return aux

atmos = np.loadtxt("falc.dat", skiprows=1)
atmos = np.array(atmos.T, dtype=np.float64)
atmos = spinor2multi(atmos)

# path to the *.input files
cwd = "."

# mu angle for which to compute the spectrum
mu = 1.0

# type of atmosphere stratification:
# 0 -- optical depth @ 500nm
# 1 -- mass density [cm2/g]
# 2 -- height [km]
atm_scale = int(0)

# wavelength samples for which to compute the spectrum (in nm)
wave_air = np.linspace(630.07, 630.3, num=251)
wave = air_to_vacuum(wave_air)

spec = pyrh.compute1d(cwd, mu, atm_scale, atmos, wave)
spec = np.array(spec, dtype=np.float64)
spec[:4] /= spec[0,0]

from astropy.io import fits
atlas = fits.open("/home/dusan/Documents/solar_atlases/Hamburg_atlas/FTS_Hamburg_center_normed_3290_8000.fits")[0].data
atlas_wave = atlas[0] / 10 # [A --> nm]
idl = np.argmin(np.abs(atlas_wave - 630.1))
atlas_int_abs = atlas[1]/atlas[1,idl] # []

plt.plot(wave_air, spec[0], color=f"C1")
plt.plot(atlas_wave, atlas_int_abs)

plt.xlabel("Wavelength in vacuum [nm]")
plt.ylabel(r"Intensity [10$^{-8}$ W/Hz/srad/m$^2$]")
plt.xlim([wave_air[0], wave_air[-1]])

plt.show()