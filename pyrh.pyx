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

	# cdef int argc
	# cdef char* argv[10]

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

	# @cython.boundscheck(False)
	# @cython.wraparound(False)
	# cpdef set_wavetable(self, cnp.ndarray[double, ndim=1, mode="c"] wave):
	# 	self.wavetable = wave
	# 	cdef Py_ssize_t i
	# 	cdef double sigma_sq
	# 	cdef double fact
	# 	self.Nwave = self.wavetable.size
	# 	for i in range(self.Nwave):
	# 		if self.wavetable[i]>199.9352:
	# 			sigma_sq = (1.0e7/self.wavetable[i])*(1.0e7/self.wavetable[i])
	# 			fact = 1.0000834213 + 2.406030e6/(1.3e10 - sigma_sq) + 1.5997e4/(3.89e9 - sigma_sq)
	# 			self.wavetable[i] = self.wavetable[i] * fact

	cpdef dummy(self):
		rh.dummy()

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef hse(self,
			  int atm_scale,
			  cnp.ndarray[double, ndim=1, mode="c"] scale,
			  cnp.ndarray[double, ndim=1, mode="c"] temp,
			  cnp.ndarray[double, ndim=1, mode="c"] ne,
			  cnp.ndarray[double, ndim=1, mode="c"] vz,
			  cnp.ndarray[double, ndim=1, mode="c"] vmic,
			  cnp.ndarray[double, ndim=1, mode="c"] mag,
			  cnp.ndarray[double, ndim=1, mode="c"] gamma,
			  cnp.ndarray[double, ndim=1, mode="c"] chi,
			  cnp.ndarray[double, ndim=2, mode="c"] nH,
			  int do_fudge,
			  cnp.ndarray[double, ndim=1, mode="c"] fudge_lam,
			  cnp.ndarray[double, ndim=2, mode="c"] fudge):
		cdef int Ndep = scale.shape[0]
		cdef int fudge_num = fudge_lam.shape[0]

		myPops = rh.hse(Ndep,
					 &scale[0], &temp[0], &ne[0], &vz[0], &vmic[0],
					 &mag[0], &gamma[0], &chi[0],
					 &nH[0,0], atm_scale,
					 do_fudge, fudge_num, &fudge_lam[0], &fudge[0,0])

		ne = convert_1d(myPops.ne, Ndep)
		nH = convert_2d(myPops.nH, 6, Ndep)

		return ne, nH

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef compute1d(self,
				int atm_scale,
				cnp.ndarray[double, ndim=1, mode="c"] scale,
				cnp.ndarray[double, ndim=1, mode="c"] temp,
				cnp.ndarray[double, ndim=1, mode="c"] ne,
				cnp.ndarray[double, ndim=1, mode="c"] vz,
				cnp.ndarray[double, ndim=1, mode="c"] vmic,
				cnp.ndarray[double, ndim=1, mode="c"] mag,
				cnp.ndarray[double, ndim=1, mode="c"] gamma,
				cnp.ndarray[double, ndim=1, mode="c"] chi,
				cnp.ndarray[double, ndim=2, mode="c"] nH,
				do_fudge,
				cnp.ndarray[double, ndim=1, mode="c"] fudge_lam,
				cnp.ndarray[double, ndim=2, mode="c"] fudge,
				cnp.ndarray[int, ndim=1, mode="c"] loggf_ids,
				cnp.ndarray[double, ndim=1, mode="c"] loggf_values,
				cnp.ndarray[int, ndim=1, mode="c"] lam_ids,
				cnp.ndarray[double, ndim=1, mode="c"] lam_values):
		cdef int Ndep = scale.shape[0]
		cdef int fudge_num = fudge_lam.shape[0]
		
		cdef int Nloggf = loggf_ids.shape[0]
		if (Nloggf!=loggf_values.shape[0]):
			print("\n  pyrh: Different number of loggf_ids and loggf_values.\n")
			sys.exit()

		cdef int Nlam = lam_ids.shape[0]
		if (Nlam!=lam_values.shape[0]):
			print("\n  pyrh: Different number of lam_ids and lam_values.\n")
			sys.exit()

		spec = rh.rhf1d(Ndep,
				 &scale[0], &temp[0], &ne[0], &vz[0], &vmic[0],
				 &mag[0], &gamma[0], &chi[0],
				 &nH[0,0], atm_scale,
				 do_fudge, fudge_num, &fudge_lam[0], &fudge[0,0],
				 Nloggf, &loggf_ids[0], &loggf_values[0],
				 Nlam, &lam_ids[0], &lam_values[0])
				 # &self.wavetable[0], self.Nwave)

		lam = convert_1d(spec.lam, spec.nlw)
		sI = convert_1d(spec.sI, spec.nlw)
		sQ, sU, sV = None, None, None
		if spec.stokes:
			sQ = convert_1d(spec.sQ, spec.nlw)
			sU = convert_1d(spec.sU, spec.nlw)
			sV = convert_1d(spec.sV, spec.nlw)
		J = convert_2d(spec.J, spec.nlw, spec.Nrays)

		Nlam = len(lam)
		return Spectrum(spec.nlw-1, lam[:Nlam-1], sI[:Nlam-1], sQ[:Nlam-1], sU[:Nlam-1], sV[:Nlam-1], J, None, spec.stokes)

	# cpdef read_RLK_lines(self):
	# 	self.rlk_lines = rh.get_RLK_lines(self.argc, self.argv)
	# 	self.Nrlk = self.rlk_lines.Nrlk

	# def get_Nrlk(self):
	# 	return self.Nrlk