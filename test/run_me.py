import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

import pyrh
from globin import rh

data = fits.open("falc_81.fits")[0].data
data = np.array(data, dtype=np.float64)

# path to the work directory
cwd = "."
cwd = '.'
# cos(theta)
mu = 1

wavelength_air = np.linspace(401.5, 401.7, num=201)
wavelength_vacuum = rh.air_to_vacuum(wavelength_air)

stokes = pyrh.compute1d(cwd, mu, 0,
						data[0,0], wavelength_vacuum,
						0, np.zeros(1, dtype=np.float64), np.zeros((3,1), dtype=np.float64),
						np.array([], dtype=np.int32), np.array([], dtype=np.float64),
						np.array([], dtype=np.int32), np.array([], dtype=np.float64))

plt.plot(wavelength_air, stokes[0])
plt.show()