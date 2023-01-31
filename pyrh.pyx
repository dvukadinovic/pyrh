# distutils: language = c
# distutils: sources = [rh/rhf1d/pyrh_compute1dray.c]
# distutils: include_dirs = rh/

cimport rh
import numpy as np
cimport numpy as cnp
import cython
import ctypes
import sys

cnp.import_array()

# cimport tools

cdef convert_1d(double *arr, Py_ssize_t n):
	cdef Py_ssize_t i
	cdef cnp.ndarray[cnp.float64_t, ndim=1] pyarr = np.empty(n)
	for i in range(n):
		pyarr[i] = arr[i]
	return pyarr

cdef convert_2d(double **arr, Py_ssize_t nx, Py_ssize_t ny):
	cdef Py_ssize_t i
	cdef Py_ssize_t j
	cdef cnp.ndarray[cnp.float64_t, ndim=2] pyarr = np.empty((nx, ny))
	for i in range(nx):
		for j in range(ny):
			pyarr[i,j] = arr[i][j]
	return pyarr

cpdef Pystring2char(lists):
	cdef int N = len(lists)
	cdef char* argv[140]

	py_string = [item.encode("utf-8") for item in lists]
	arr = (ctypes.c_char_p * N)(*py_string)
	for i_ in range(N):
		argv[i_] = arr[i_]

	return N, argv

class Spectrum(object):
	def __init__(self, nlw=None, lam=None, sI=None, sQ=None, 
				 sU=None, sV=None, J=None, Jrh=None, stokes=True):
		self.nlw = nlw
		self.lam = lam
		self.I = sI
		self.Q = None
		self.U = None
		self.V = None
		self.J = J
		self.Jrh = Jrh
		self.stokes = False

		if stokes:
			self.stokes = True
			self.Q = sQ
			self.U = sU
			self.V = sV
		
# add OF table as input to rhf1d()

cdef class RH:

	# cdef int Nrlk
	# cdef rh.myRLK_Line rlk_lines

	# cdef Py_ssize_t Nwave
	# cdef array[double] wavetable
	# cdef double[::1] wavetable

	def __init__(self):
		pass

	# def __reduce__(self):
	# 	return (self.__class__, (None, None))

	# def __init__(self, input="keyword.input", logfile=None, quiet=True):
	# 	py_argv = "rhf1d" + " -input " + input
	# 	if logfile is not None:
	# 		py_argv += " -logfile " + logfile
	# 	if quiet:
	# 		py_argv += " -quiet"

	# 	py_list = py_argv.split(" ")
	# 	self.argc = len(py_list)
	# 	py_string = [item.encode("utf-8") for item in py_list]
	# 	arr = (ctypes.c_char_p * self.argc)(*py_string)
	# 	cdef Py_ssize_t i_
	# 	for i_ in range(self.argc):
	# 		self.argv[i_] = arr[i_]

		# self.Nrlk = 0
		# self.Nwave = 0

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef set_wavetable(self, cnp.ndarray[double, ndim=1, mode="c"] wave):
		cdef Py_ssize_t i
		cdef double sigma_sq
		cdef double fact
		cdef Py_ssize_t Nwave = wave.size
		for i in range(Nwave):
			if wave[i]>199.9352:
				sigma_sq = (1.0e7/wave[i])*(1.0e7/wave[i])
				fact = 1.0000834213 + 2.406030e6/(1.3e10 - sigma_sq) + 1.5997e4/(3.89e9 - sigma_sq)
				wave[i] = wave[i] * fact

		return wave

	# cpdef dummy(self):
	# 	rh.dummy()

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef hse(cwd,
		  int atm_scale,
		  cnp.ndarray[double, ndim=1, mode="c"] scale,
		  cnp.ndarray[double, ndim=1, mode="c"] temp,
		  double pg_top,
		  int do_fudge,
		  cnp.ndarray[double, ndim=1, mode="c"] fudge_lam,
		  cnp.ndarray[double, ndim=2, mode="c"] fudge):
	cdef int Ndep = scale.shape[0]
	cdef int fudge_num = fudge_lam.shape[0]

	cdef char* argv[140]

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	myPops = rh.hse(argv[0], Ndep, pg_top,
				 &scale[0], &temp[0],
				 atm_scale,
				 do_fudge, fudge_num, &fudge_lam[0], &fudge[0,0])

	ne = convert_1d(myPops.ne, Ndep)
	nHtot = convert_1d(myPops.nHtot, Ndep)
	nH = convert_2d(myPops.nH, 6, Ndep)
	rho = convert_1d(myPops.rho, Ndep)
	pg = convert_1d(myPops.pg, Ndep)

	return ne, nH, nHtot, rho, pg

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef get_ne_from_nH(cwd,
					 int atm_scale,
					 cnp.ndarray[double, ndim=2, mode="c"] atmosphere,
					 cnp.ndarray[double, ndim=1, mode="c"] nH):
	cdef int Ndep = atmosphere.shape[1]

	cdef char* argv[140]

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	rh.get_ne_from_nH(argv[0], atm_scale, Ndep, 
			  &atmosphere[0,0], &atmosphere[1,0], 
			  &nH[0], &atmosphere[2,0])

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef get_tau(cwd,
			  double mu,
			  int atm_scale,
			  cnp.ndarray[double, ndim=2, mode="c"] atmosphere,
			  double lam_ref):
	cdef int Ndep = atmosphere.shape[1]

	cdef char* argv[140]

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	cdef cnp.ndarray[double, ndim=1, mode="c"] tau = np.ones(Ndep)

	rh.get_tau(argv[0], mu, Ndep, &tau[0],
			 &atmosphere[0,0], &atmosphere[1,0], 
			 &atmosphere[2,0], &atmosphere[3,0], &atmosphere[4,0],
			 &atmosphere[5,0], &atmosphere[6,0], &atmosphere[7,0],
			 &atmosphere[8,0], atm_scale,
			 lam_ref)

	return tau

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef compute1d(cwd,
				double mu,
				int atm_scale,
				cnp.ndarray[double, ndim=2, mode="c"] atmosphere,
				cnp.ndarray[double, ndim=1, mode="c"] wave,
				int do_fudge,
				cnp.ndarray[double, ndim=1, mode="c"] fudge_lam,
				cnp.ndarray[double, ndim=2, mode="c"] fudge,
				cnp.ndarray[int, ndim=1, mode="c"] loggf_ids,
				cnp.ndarray[double, ndim=1, mode="c"] loggf_values,
				cnp.ndarray[int, ndim=1, mode="c"] lam_ids,
				cnp.ndarray[double, ndim=1, mode="c"] lam_values):
	cdef int Ndep = atmosphere.shape[1]
	cdef int fudge_num = fudge_lam.shape[0]
	cdef int Nwave = wave.shape[0]
	
	cdef int Nloggf = loggf_ids.shape[0]
	if (Nloggf!=loggf_values.shape[0]):
		print("\n  pyrh: Different number of loggf_ids and loggf_values.\n")
		sys.exit()

	cdef int Nlam = lam_ids.shape[0]
	if (Nlam!=lam_values.shape[0]):
		print("\n  pyrh: Different number of lam_ids and lam_values.\n")
		sys.exit()

	cdef char* argv[140]

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	spec = rh.rhf1d(argv[0], mu, Ndep,
			 &atmosphere[0,0], &atmosphere[1,0], 
			 &atmosphere[2,0], &atmosphere[3,0], &atmosphere[4,0],
			 &atmosphere[5,0], &atmosphere[6,0], &atmosphere[7,0],
			 &atmosphere[8,0], atm_scale,
			 Nwave, &wave[0],
			 do_fudge, fudge_num, &fudge_lam[0], &fudge[0,0],
			 Nloggf, &loggf_ids[0], &loggf_values[0],
			 Nlam, &lam_ids[0], &lam_values[0],
			 0, argv[0])
			 # &self.wavetable[0], self.Nwave)

	# spec.nlw -= 1
	lam = convert_1d(spec.lam, spec.nlw)
	sI = convert_1d(spec.sI, spec.nlw)
	sQ, sU, sV = None, None, None
	if spec.stokes:
		sQ = convert_1d(spec.sQ, spec.nlw)
		sU = convert_1d(spec.sU, spec.nlw)
		sV = convert_1d(spec.sV, spec.nlw)

	# J = convert_2d(spec.J, spec.nlw, spec.Nrays)

	# Nlam = len(lam)
	return sI, sQ, sU, sV
	# return Spectrum(spec.nlw, lam, sI, sQ, sU, sV, None, None, spec.stokes)

	# cpdef read_RLK_lines(self):
	# 	self.rlk_lines = rh.get_RLK_lines(self.argc, self.argv)
	# 	self.Nrlk = self.rlk_lines.Nrlk

	# def get_Nrlk(self):
	# 	return self.Nrlk

