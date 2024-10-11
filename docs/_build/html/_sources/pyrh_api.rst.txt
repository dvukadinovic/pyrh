pyrh API
========

.. py:function:: pyrh.hse(cwd, atm_scale, scale, temp, pg_top=0.1, fudge_wave=None, fudge_value=None, full_output=False)

	Set the atmosphere in hydrostatic equilibrium (HSE) by computing the electron and total hydrogen densities from given temperature, assuming an ideal gas law and LTE populations of species in the atmosphere.

	:param cwd: Relative path to directory containing RH input files.
	:type cwd: str

	:param atm_scale: Type of an atmosphere scale. ``atm_scale=0`` for the optical depth, ``atm_scale=1`` for the column mass density and ``atm_scale=2`` for the height.
	:type atm_scale: int

	:param scale: Array containing the atmosphere scale.
	:type scale: C contiguouse float64 numpy.ndarray

	:param temp: Array containing a temperature in K for each depth point specified in ``scale``.
	:type temp: C contiguouse float64 numpy.ndarray

	:param pg_top: Gas pressure at the top of the atmosphere in Pa. It is used to start the HSE iterative solution. By default ``pg_top=0.1``.
	:type pg_top: default, double

	:param fudge_wave: Wavelength grid for the fudge factors used to alter the continuum opacity. By default ``fudge_wave=None``.
	:type fudge_wave: default, C contiguouse float64 numpy.ndarray

	:param fudge_value: The opacity fudge coefficients for each wavelength given in ``fudge_wave``. This is 2D array with three columns, each giving the fudge coefficients for H-, scattering and metals continuum opacity sources, respectively. By default ``fudge_value=None``.
	:type fudge_value: default, C contiguouse float64 numpy.ndarray

	:param full_output: Flag for the amount of output provided. If ``True``, the method return electron density, total hydrogen density, mass density and gas pressure, respectively. By default ``full_output=False`` and the method returns only electron density and total hydrogen density.
	:type full_output: default, bool

	:return: Electron and total hydrogen density. If ``full_output=True``, then it returns electron density, total hydrogen density, mass density and gas pressure.

.. py:function:: pyrh.get_scales(cwd, atm_scale, scale, atmosphere, lam_ref)

	Compute the atmospheric scales: optical depth, height and column mass density at a given reference wavelength. From the provided scale, method computes opacity and converts it to the other two scales.

	:param cwd: Relative path to directory containing RH input files.
	:type cwd: str

	:param atm_scale: Type of an atmosphere scale. ``atm_scale=0`` for the optical depth, ``atm_scale=1`` for the column mass density and ``atm_scale=2`` for the height.
	:type atm_scale: int

	:param scale: Array containing the atmosphere scale.
	:type scale: C contiguouse float64 numpy.ndarray

	:param atmosphere: MULTI type atmosphere (check ... for more details).
	:type atmosphere: C contiguouse float64 numpy.ndarray

	:param lam_ref: Reference wavelength in vacuum in nanometer units for which we give scales (optical depth and column mass) or at which we want to compute them.
	:type lam_ref: double

	:return: The scales in this order: optical depth, height in meters and column mass in kg/m2.
	:return type: C contiguouse float64 numpy.ndarray

.. py:function:: pyrh.get_ne_from_nH(cwd, atm_scale, scale, temperature, nH)

	Compute the LTE electron density from given temperature and total hydrogen density. 

	:param cwd: Relative path to directory containing RH input files.
	:type cwd: str

	:param atm_scale: Type of an atmosphere scale. ``atm_scale=0`` for the optical depth, ``atm_scale=1`` for the column mass density and ``atm_scale=2`` for the height.
	:type atm_scale: int

	:param scale: Array containing the atmosphere scale.
	:type scale: C contiguouse float64 numpy.ndarray

	:param temp: Array containing a temperature in K for each depth point specified in ``scale``.
	:type temp: C contiguouse float64 numpy.ndarray

	:param nH: The total hydrogen density in 1/cm3.
	:type nH: C contiguouse float64 numpy.ndarray

	:return: The electron density in 1/cm3.
	:return type: C contiguouse float64 numpy.ndarray

.. py:function:: pyrh.compute1d(cwd, mu, atm_scale, atmosphere, wave, loggf_ids=None, loggf_values=None, lam_ids=None, lam_values=None, fudge_wave=None, fudge_value=None, get_atomic_rfs=False)

	Synthesise a spectrum and return the Stokes vector.

	:param cwd: Relative path to directory containing RH input files.
	:type cwd: str

	:param mu: Angle for which we are computing the spectrum. 
	:type mu: double

	:param atm_scale: Type of an atmosphere scale. ``atm_scale=0`` for the optical depth, ``atm_scale=1`` for the column mass density and ``atm_scale=2`` for the height.
	:type atm_scale: int

	:param atmosphere: MULTI type atmosphere (check ... for more details).
	:type atmosphere: C contiguouse float64 numpy.ndarray

	:param wave: Wavelength in vacuum in nanometer units for which to synthesise a spectrum.
	:type wave: C contiguouse float64 numpy.ndarray

	:param loggf_ids: Spectral line number from a Kurucz line list for which we provide log(gf) value different from the one found in a line list. Default ``None``.
	:type loggf_ids: C contiguouse float64 numpy.ndarray

	:param loggf_values: log(gf) values for each line from ``loggf_ids``. Default ``None``.
	:type loggf_values: C contiguouse float64 numpy.ndarray

	:param lam_ids: Spectral line number from a Kurucz line list for which we alter the central wavelength. Default ``None``.
	:type lam_ids: C contiguouse float64 numpy.ndarray

	:param lam_values: The wavelength shift in mA for all lines given in ``lam_ids``.
	:type lam_values: C contiguouse float64 numpy.ndarray

	:param fudge_wave: Wavelength grid for the fudge factors used to alter the continuum opacity. By default ``fudge_wave=None``.
	:type fudge_wave: default, C contiguouse float64 numpy.ndarray

	:param fudge_value: The opacity fudge coefficients for each wavelength given in ``fudge_wave``. This is 2D array with three columns, each giving the fudge coefficients for H-, scattering and metals continuum opacity sources, respectively. By default ``fudge_value=None``.
	:type fudge_value: default, C contiguouse float64 numpy.ndarray

	:param get_atomic_rfs: Flag if we want to compute the analytical response functions to log(gf) for lines given in ``loggf_ids``. Default ``False``.
	:type get_atomic_rfs: bool

	:return: The Stokes vector and wavelength in vacuum in nanometers.