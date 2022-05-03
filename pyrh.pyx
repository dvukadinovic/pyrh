# distutils: language = c
# distutils: sources = [rh/rhf1d/pyrh_compute1dray.c]
# distutils: include_dirs = rh/

cimport rh
import ctypes
import numpy as np
cimport numpy as cnp
import cython
from libc.stdlib cimport malloc, free

# cimport tools

cdef convert_1d(double *arr, int n):
	cdef int i
	cdef cnp.ndarray[cnp.float64_t, ndim=1] pyarr = np.empty(n)
	for i in range(n):
		pyarr[i] = arr[i]
	return pyarr

cdef convert_2d(double **arr, int nx, int ny):
	cdef int i
	cdef int j
	cdef cnp.ndarray[cnp.float64_t, ndim=2] pyarr = np.empty((nx, ny))
	for i in range(nx):
		for j in range(ny):
			pyarr[i,j] = arr[i][j]
	return pyarr

class Spectrum(object):
	def __init__(self, nlw, lam, sI, sQ, sU, sV, J, Jrh, stokes):
		self.nlw = nlw
		self.lam = lam
		self.I = sI
		self.Q = sQ
		self.U = sU
		self.V = sV
		self.J = J
		self.Jrh = Jrh

		if stokes:
			self.stokes = True
		else:
			self.stokes = False

cdef class RH:
	cdef rh.mySpectrum spec

	cdef int Nrlk
	cdef rh.myRLK_Line rlk_lines

	cdef int argc
	cdef char* argv[10]

	cdef int Nwave
	cdef double* wavetable_ptr

	def __init__(self, input="keyword.input", logfile=None, quiet=True):
		py_argv = "rhf1d" + " -input " + input
		if logfile is not None:
			py_argv += " -logfile " + logfile
		if quiet:
			py_argv += " -quiet"

		py_list = py_argv.split(" ")
		self.argc = len(py_list)
		py_string = [item.encode("utf-8") for item in py_list]
		arr = (ctypes.c_char_p * self.argc)(*py_string)
		cdef int i_
		for i_ in range(self.argc):
			self.argv[i_] = arr[i_]

		self.Nrlk = 0
		self.Nwave = 0

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef set_wavetable(self, cnp.ndarray[double, ndim=1, mode="c"] wavetable):
		cdef int i
		cdef double sigma_sq
		cdef double fact
		self.Nwave = wavetable.size
		# self.wavetable = (double *)malloc(self.Nwave * cython.sizeof(double))
		for i in range(self.Nwave):
			if wavetable[i]>199.9352:
				sigma_sq = (1.0e7/wavetable[i])*(1.0e7/wavetable[i])
				fact = 1.0000834213 + 2.406030e6/(1.3e10 - sigma_sq) + 1.5997e4/(3.89e9 - sigma_sq)
				wavetable[i] = wavetable[i] * fact
		self.wavetable_ptr = &wavetable[0]

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef compute1d(self,
				cnp.ndarray[double, ndim=1, mode="c"] scale,
				cnp.ndarray[double, ndim=1, mode="c"] temp,
				cnp.ndarray[double, ndim=1, mode="c"] ne,
				cnp.ndarray[double, ndim=1, mode="c"] vz,
				cnp.ndarray[double, ndim=1, mode="c"] vmic,
				cnp.ndarray[double, ndim=1, mode="c"] mag,
				cnp.ndarray[double, ndim=1, mode="c"] gamma,
				cnp.ndarray[double, ndim=1, mode="c"] chi,
				cnp.ndarray[double, ndim=2, mode="c"] nH,
				int atm_scale):
		cdef int Ndep = scale.shape[0]

		self.spec = rh.rhf1d(self.argc, self.argv, Ndep,
				 &scale[0], &temp[0], &ne[0], &vz[0], &vmic[0],
				 &mag[0], &gamma[0], &chi[0],
				 &nH[0,0], atm_scale, &self.rlk_lines,
				 self.wavetable_ptr, self.Nwave)

		lam = convert_1d(self.spec.lam, self.spec.nlw)
		sI = convert_1d(self.spec.sI, self.spec.nlw)
		sQ, sU, sV = None, None, None
		if self.spec.stokes:
			sQ = convert_1d(self.spec.sQ, self.spec.nlw)
			sU = convert_1d(self.spec.sU, self.spec.nlw)
			sV = convert_1d(self.spec.sV, self.spec.nlw)
		J = convert_2d(self.spec.J, self.spec.nlw, self.spec.Nrays)

		#lam_ = np.ctypeslib.as_array(self.spec.lam, shape=(self.spec.nlw,))

		return Spectrum(self.spec.nlw, lam, sI, sQ, sU, sV, J, None, self.spec.stokes)

	cpdef read_RLK_lines(self):
		self.rlk_lines = rh.get_RLK_lines(self.argc, self.argv)
		self.Nrlk = self.rlk_lines.Nrlk

	def get_Nrlk(self):
		return self.Nrlk