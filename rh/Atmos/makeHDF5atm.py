"""

Script for converting MURaM atmosphere from fits to HDF5 format neccessary
for running RH15 code.

"""
from helita.sim.rh15d import make_xarray_atmos
from astropy.io import fits
import numpy as np
from scipy.constants import k
import sys

def hydrogen_lvl_pops(nH, temp, nlvl=6):
	"""
	Redistribute hydrogen atoms through given number of levels ('nlvl').

	Parameters:
	---------------
	nH : ndarray
		total hydrogen population through atmosphere (1D or 3D array).
	temp : ndarray
		temperature stratification in atmosphere (same dimenstion as nH).
	nlvl : int (optional)
		number of levels in hydrogen atom for which to calculate populations.

	Return:
	---------------
	popos : ndarray
		populations of hydrogen levels. Dimension is (nvlv, shape(nH)).
	"""
	pops = np.zeros((nlvl, *nH.shape))

	for lvl in range(nlvl):
		e_lvl = 13.6*(1-1/(lvl+1)**2)
		pops[lvl] = nH/2 * 2*(lvl+1)**2 * np.exp(-5040/temp * e_lvl)
	return pops

# element abundances specified by RH code
A = np.loadtxt("../Atoms_example/abundance.input", usecols=(1,), unpack=True)

# input atmosphere
# np=12 x nx=288 x ny=288 x nz=100
atmos = fits.open("/data/slam/home/vukadinovic/cubes/muram/muram_50G.fits")[0].data
shape = atmos.shape

# total number of gas particles (excluding electrons)
ng = (atmos[3]-atmos[4])*0.1/k/atmos[2]
abund = np.sum(10**(A-12))
# Phi function from Saha equation for hydrogen atom
phi = 0.6665 * 1/2 * atmos[2]**(5/2) * 10**(-5040/atmos[2]*13.6)
# hydrogen level populations
nH = hydrogen_lvl_pops(ng/abund, atmos[2])

temp = atmos[2]		   # [K]
vz = atmos[9] / 1e2    # [m/s]
z = atmos[1] / 1e2     # [m]
rho = atmos[6]* 1e3    # [kg / m3]
vturb = atmos[8] / 1e5 # [km/s]
ne = atmos[4]*0.1 / k / atmos[2] # [1/m3]
Bx = atmos[7]*np.sin(atmos[10])*np.cos(atmos[11]) / 1e4 # [T]
By = atmos[7]*np.sin(atmos[10])*np.sin(atmos[11]) / 1e4 # [T]
Bz = atmos[7]*np.cos(atmos[10]) / 1e4 # [T]
x = np.arange(0,shape[1],1, dtype=np.int)
y = np.arange(0,shape[2],1, dtype=np.int)

description = "MURaM atmosphere cube of quiet Sun with 50G mean field."

# save atmospheric model to HDF5 file format
make_xarray_atmos(outfile="muram_50G.hdf5",
	T=temp, vz=vz, 
	x=x, y=y, z=z, rho=rho, ne=ne, nH=nH,
	vturb=vturb, Bx=Bx, By=By, Bz=Bz, 
	desc=description)

