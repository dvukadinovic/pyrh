import numpy as np
import matplotlib.pyplot as plt
import pyrh
from scipy.interpolate import interp1d
from astropy.io import fits
from globin.atmos import write_multi_atmosphere
from globin.atmos import distribute_hydrogen
from globin.rh import write_wavs

def spinor2multi(atm):
	"""
	Casting from SPINOR type atmosphere structure to MULTI type structure.
	"""
	from scipy.constants import k
	new = np.empty((14, atm.shape[-1]), dtype=np.float64)

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

def vacuum_to_air(wavelength, vacuum2air_limit=200.0000):
    factor = np.ones_like(wavelength)
    wave2 = 1/wavelength**2
    factor[wavelength>vacuum2air_limit] = 1 + 2.735182e-4 + (1.314182 + 2.76249e4*wave2) * wave2

    aux = wavelength/factor

    return aux

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

# mu value for which to compute the spectrum
mu = 1.0

# type of atmosphere stratification:
# 0 -- optical depth @ 500nm
# 1 -- mass density
# 2 -- height [km]
atm_scale = int(0)

# wavelength grid in vacuum on which we compute the spectrum (in nm)
wave = np.linspace(630.2-0.15, 630.2+0.15, num=301)
wave_vacuum = air_to_vacuum(wave)

spec = pyrh.compute1d(cwd, mu, atm_scale, atmos,
					  wave_vacuum,
					  int(0), np.empty(1, dtype=np.float64), np.empty((3,1), dtype=np.float64),
					  np.array([], dtype=np.int32), np.array([], dtype=np.float64),
					  np.array([], dtype=np.int32), np.array([], dtype=np.float64),
					  False)
spec = np.array(spec, dtype=np.float64)
# np.savetxt("hinode_falc_pyrh.spec", spec.T, fmt="%7.6e")

#--- create figure
fig, ax = plt.subplots(nrows=2, ncols=1, 
	sharex=True, height_ratios=[4,1])
plt.subplots_adjust(hspace=0)

#--- spectra
spec[0] /= spec[0,0]
ax[0].plot(wave, spec[0], label="synth", c="k")

hinode_pyrh = np.loadtxt("hinode_falc_pyrh.spec").T
hinode_pyrh[0] /= hinode_pyrh[0,0]
ax[0].plot(wave, hinode_pyrh[0], label="pyrh benchamark", c="tab:red")

hinode_rh = np.loadtxt("hinode_falc_rh.spec").T
hinode_rh[0] /= hinode_rh[0,0]
ax[0].plot(wave, hinode_rh[0], label="rh benchamark", c="tab:green")

# spinor = fits.open("inverted_profs.1.fits")[0].data[0,0]
# spinor /= spinor[0,1]
# ax[0].plot(wave, spinor[0], label="SPINOR benchamark", c="tab:blue")

#--- differences
ax[1].plot(wave, spec[0] - hinode_pyrh[0], c="tab:red")
ax[1].plot(wave, spec[0] - hinode_rh[0], c="tab:green")
# ax[1].plot(wave, spec[0] - spinor[0], c="tab:blue")

ax[0].set_ylabel(r"Normalized intensity")
ax[1].set_xlabel("Wavelength in air [nm]")
ax[1].set_xlim([wave[0], wave[-1]])

ax[0].legend(frameon=True, 
	ncols=4, 
	bbox_to_anchor=(0, 1.01, 1, 0.2),)

plt.show()