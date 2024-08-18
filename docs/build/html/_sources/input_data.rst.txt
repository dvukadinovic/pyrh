Input parameters
================

Along the necesary input files used by RH, we are going to describe here the structure of an atmospheric model that is passed to *pyrh* and the wavelength vector for which we are synthesising a spectrum. Other input parameters of less imnportance are described along with each *pyrh* method.

Atmosphere
^^^^^^^^^^

The main ingredient to compute spectra is an atmospheric model. The structure of the atmospheric model used for spectral synthesis using *pyrh* is the MULTI type atmosphere which has the following order of atmospheric parameters:

.. table:: Order of atmospheric parameters.

	========================== =======
	parameter                  unit
	========================== =======
	scale                      none, km or cm2/g
	temperature                K 
	electron density           1/cm3
	line-of-sight velocity     km/s
	micro-turbulent velocity   km/s
	magnetic field strength    G
	magnetic field inclination deg
	magnetic field azimuth     deg
	total hydrogen density     1/cm3
	========================== =======

The scale of an atmosphere can be given in height (in km), optical depth (assumed to be at 500 nm) or column mass density (cm2/g). If the scale is given in the column mass density, total hydrogen density will be computed internally by RH and the input will be disregarded. 

The type of the scale of an atmosphere needs to be passed as an argument to any method from *pyrh*. We have ``atm_scale=0`` for the optical depth scale, ``atm_scale=1`` for the column mass density and ``atm_scale=2`` for the height.

The top of the atmosphere is assumed to be the first row, while the bottom is the last row.

Wavelength vector
^^^^^^^^^^^^^^^^^

To compute a spectrum we also need to provide an array containing wavelengths in nanometers. These wavelengths are assumed to be given in vacuum. 

Internally, RH will add a new wavelength point to array that refers to the reference wavelength at which the optical depth scale is calculated. By default this is 500 nm unless otherwise specified by the user in ``keyword.input`` file.

The output Stokes vector is given in absolute units W/Hz/srad/m2.