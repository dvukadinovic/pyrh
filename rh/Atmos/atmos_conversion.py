import numpy as np
from scipy.constants import k
from astropy.io import fits
from scipy.interpolate import interp1d

def hydrogen_lvl_pops(temp, pg, pe, nlvl=6):
	"""
	Redistribute hydrogen atoms in first nlvl-1. Last column
	is reserved for proton number density.

	Parameters:
	---------------
	temp : ndarray
		temperature stratification in atmosphere (same dimenstion as nH).
	pg : ndarray
		gas pressure (1D or 3D array).
	pe : ndarray
		electron pressure (1D or 3D array).
	nlvl : int (optional)
		number of levels in hydrogen atom for which to calculate populations.
		last index stores proton numbers.

	Return:
	---------------
	popos : ndarray
		populations of hydrogen levels + protons. Dimension is (nvlv, shape(nH)).
	"""
	import matplotlib.pyplot as plt

	nH = (pg-pe)/10 / k/temp / np.sum(10**(abundance-12)) / 1e6
	nH0 = nH / (1 + saha_phi(temp)/pe)
	nprot = nH - nH0
	
	tt = np.linspace(3000,10000,num=8)
	u0_coeffs = [2.00000e+00, 2.00000e+00, 2.00000e+00, 2.00000e+00, 2.00000e+00, 2.00001e+00, 2.00003e+00, 2.00015e+00]
	U0 = interp1d(tt, u0_coeffs)(temp)

	pops = np.zeros((nlvl, *nH.shape))
	pops[-1] = nprot
	for lvl in range(nlvl-1):
		e_lvl = 13.59844*(1-1/(lvl+1)**2)
		g = 2*(lvl+1)**2
		pops[lvl] = nH0/U0 * g * np.exp(-5040/temp * e_lvl)

	return pops

def saha_phi(
	temp, 
	u0_coeffs=[2.00000e+00, 2.00000e+00, 2.00000e+00, 2.00000e+00, 2.00000e+00, 2.00001e+00, 2.00003e+00, 2.00015e+00], 
	u1_coeffs=[1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00],
	Ej=13.59844):
	"""
	Calculate Phi(T) function for Saha's equation in form:

	n+/n0 = Phi(T)/Pe

	All units are in cgs system.

	We interpolate for partition function coefficient. Temperature sale on which
	values are expected are from 3000K to 10000K with step of 1000K.

	Parameters:
	---------------
	temp : ndarray
		temperature for which to calculate Phi(T) function
	u0_coeffs : list
		partition function of lower ionization stage from 3000K to 10000K.
	u1_coeffs : list
		partition function of higher ionization stage
	Ej : float
		ionization energy of state in [eV]

	Return:
	---------------
	Phi(T) : ndarray
		value of Phi(T) function at every temperature
	"""
	u1 = interp1d(np.linspace(3000,10000,num=8),u1_coeffs)(temp)
	u0 = interp1d(np.linspace(3000,10000,num=8),u0_coeffs)(temp)
	return 0.6665 * u1/u0 * temp**(5/2) * 10**(-5040/temp*Ej)

def interpolate_cube(cube, x):
	shape = cube.shape
	new_cube = np.zeros((*shape[:3],len(x)))

	for idx in range(shape[0]):
		for idy in range(shape[1]):
			new_cube[idx,idy,0,:] = x
			for pID in range(1,shape[2]):
				fun = interp1d(cube[idx,idy,0], cube[idx,idy,pID])
				new_cube[idx,idy,pID] = fun(x)

	return new_cube

nlvl = 6
abundance = np.loadtxt("../Atoms/abundance.input", usecols=(1,), unpack=True)

logt, temp, pg, pe, vturb, vz = np.loadtxt("falc.dat", 
							skiprows=1, unpack=True, usecols=(0,2,3,4,8,9))
ne = pe/10/k/temp / 1e6 # [m-3]
vturb /= 1e5 # [km/s]

n_pops = hydrogen_lvl_pops(temp, pg, pe, nlvl)

#--- make single atmos output model
out = open("FALC_tau_scale.atmos", "w")

out.write("* Model file created on 23/10/2020 from 'atmos_conversion.py'\n")
out.write("*\n")
out.write("  FALC_tau_scale\n")
out.write("  Tau scale\n")
out.write("*\n")
out.write("* log(g) [cm s^-2]\n")
out.write("  4.44\n")
out.write("*\n")
out.write("* Ndep\n")
out.write(f"  {len(temp)}\n")
out.write("*\n")
out.write("* lot tau    Temperature        Ne         V              Vturb\n")

for i_ in range(len(temp)):
	out.write("  {:+5.4f}   {:6.2f}   {:5.4e}   {:5.4e}   {:5.4e}\n".format(logt[i_], temp[i_], ne[i_], vz[i_], vturb[i_]))

out.write("*\n")
out.write("* Hydrogen populations\n")
out.write("*  nh(1)       nh(2)       nh(3)       nh(4)       nh(5)       nh(6)\n")

for i_ in range(len(temp)):
	out.write("   {:5.4e}   {:5.4e}   {:5.4e}   {:5.4e}   {:5.4e}   {:5.4e}\n".format(n_pops[0,i_], n_pops[1,i_], n_pops[2,i_], n_pops[3,i_], n_pops[4,i_], n_pops[5,i_]))

out.close()

#--- make fits file format of atmosphere
nx, ny = 2,3
atmos_cube = np.zeros((nx,ny,14,len(temp)))

# inv_vz = 0.00064134 # [km/s]
vz = 0 # [km/s]

# standard physical parameters
atmos_cube[:,:,0,:] = logt
atmos_cube[:,:,1,:] = temp
atmos_cube[:,:,2,:] = ne
atmos_cube[:,:,3,:] = vz
atmos_cube[:,:,4,:] = vturb
# magnetic field strength
atmos_cube[:,:,5,:] = 250 / 1e4 # G --> T
# magnetic field inclination
atmos_cube[:,:,6,:] = np.deg2rad(45)
# magnetic field azimuth
atmos_cube[:,:,7,:] = np.deg2rad(45)

# hydrogen populations by levels
for i_ in range(nlvl):
	atmos_cube[:,:,8+i_,:] = n_pops[i_]

# logtau = np.linspace(-6,1,num=71)
# cube = interpolate_cube(atmos_cube, logtau)

print(atmos_cube.shape)

fits.writeto("falc_1x1.fits", atmos_cube, overwrite=True)