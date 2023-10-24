import numpy as np
import matplotlib.pyplot as plt
import pyrh
from scipy.interpolate import interp1d

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

# mu value for which to compute the spectrum
mu = 1.0

# type of atmosphere stratification:
# 0 -- optical depth @ 500nm
# 1 -- mass density
# 2 -- height [km]
atm_scale = int(0)

# wavelength grid in vacuum on which we compute the spectrum (in nm)
wave = np.linspace(630.0, 630.35, num=251) + 0.2

spec = pyrh.compute1d(cwd, mu, atm_scale, atmos,
					  wave,
					  int(0), np.empty(1, dtype=np.float64), np.empty((3,1), dtype=np.float64),
					  np.array([], dtype=np.int32), np.array([], dtype=np.float64),
					  np.array([], dtype=np.int32), np.array([], dtype=np.float64),
					  False)
spec = np.array(spec, dtype=np.float64)

fig, ax = plt.subplots(nrows=2, ncols=1, 
	sharex=True, height_ratios=[4,1])
plt.subplots_adjust(hspace=0)

#--- spectra

ax[0].plot(wave, spec[0], 
	label="synth", 
	c="k")

hinode_pyrh = np.loadtxt("hinode_falc_pyrh.spec").T
ax[0].plot(hinode_pyrh[-1], hinode_pyrh[0], 
	label="pyrh benchamark",
	c="tab:red")

# hinode_pyrh = np.loadtxt("hinode_falc_spinor.spec").T
# ax[0].plot(hinode_pyrh[-1], hinode_pyrh[0], label="SPINOR benchamark")

#--- differences
x = hinode_pyrh[-1]
y = hinode_pyrh[0]
if not np.allclose(wave, hinode_pyrh[-1], atol=1e-6, equal_nan=True):
	x = wave
	y = interp1d(hinode_pyrh[-1], y, kind="cubic")(x)
ax[1].plot(x, spec[0] - y,
	c="tab:red")

ax[0].set_ylabel(r"Intensity [W Hz$^{-1}$ sr$^{-1}$]")
ax[1].set_xlabel("Wavelength in vacuum [nm]")
ax[1].set_xlim([wave[0], wave[-1]])

ax[0].legend(frameon=True, loc="best")

plt.show()