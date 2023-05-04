Quick start examples
================================

By default `pyrh' is created as a Python gateway to the C code of RH for spectrum synthesis. We will now show how can we compute some spectra using `pyrh`. THe user should be familiar with the RH code in general since the `pyrh` does not provide any configuration settings. The basic usage is done through ``*.input`` files of RH and user should use these as he would use RH to synthesize spectra.

Open the ``synth.py`` folder from root directory of `pyrh` module. 

We first read in the atmosphere. Here we use FAL C atmoshere model:

.. code-block:: Python
	:linenos:

	import numpy as np
	import matplotlib.pyplot as plt
	import pyrh

	def spinor2multi(atm):
		"""
		Casting from SPINOR type atmosphere structure to MULTI type structure.
		"""
		from scipy.constants import k
		new = np.empty(atm.shape, dtype=np.float64)

		new[0] = atm[0]
		new[1] = atm[2]
		# electron number dnesity
		new[2] = atm[4]/10/k/atm[2] / 1e6 # [1/cm3]
		new[3] = atm[9]/1e5
		new[4] = atm[8]/1e5
		# magnetic field vector: strength, inclination and azimuth
		new[5] = atm[7]
		new[6] = atm[-2]
		new[7] = atm[-1]
		# total hydrogen number density
		new[8] = (atm[3] - atm[4])/10/k/atm[2] / 1e6 / 1.26 # [1/cm3]

		return new

	atmos = np.loadtxt("falc.dat", skiprows=1)
	atmos = np.array(atmos.T, dtype=np.float64)
	atmos= spinor2multi(atmos)

We here had to convert the atmosphere from SPINOR like format found in the `falc.dat` file to the MULTI type atmosphere which is used by RH.

Then, we define couple of variables that are used as configurations for RH:

.. code-block:: Python
	:linenos:

	# path to the *.input files
	cwd = "."

	# mu angle for which to compute the spectrum
	mu = 1.0

	# type of atmosphere stratification:
	# 0 -- optical depth @ 500nm
	# 1 -- mass density
	# 2 -- height [km]
	atm_scale = int(0)

Further, we set up the wavelength grid for which we are synthesizing the spectrum. We will synthesize Fe I 6301 and 6302 line pair, so:

.. code-block:: Python
	:linenos:

	# wavelength samples for which to compute the spectrum (in nm) in vacuum
	wave = np.linspace(630.0, 630.35, num=251) + 0.2

Now, we call the ``pyrh.compute1d`` which computes the spectrum for given input parameters.

.. code-block:: Python
	:linenos:

	spec = pyrh.compute1d(cwd, mu, atm_scale, atmos,
		wave,
		int(0), np.empty(1, dtype=np.float64), np.empty((3,1), dtype=np.float64),
		np.array([], dtype=np.int32), np.array([], dtype=np.float64),
		np.array([], dtype=np.int32), np.array([], dtype=np.float64))
	
And the final product is:

.. code-block:: Python
	:linenos:

	spec = np.array(spec, dtype=np.float64)

	plt.plot(wave, spec[0])
	plt.show()

Voila!