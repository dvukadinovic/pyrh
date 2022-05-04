# distutils: language = c
# distutils: sources = [rh/rhf1d/pyrh_compute1dray.c]
# distutils: include_dirs = rh/

cimport rh
import ctypes
import numpy as np
cimport numpy as cnp
import cython
from libc.stdlib cimport malloc, free

from cpython cimport PyObject, Py_INCREF

cnp.import_array()

# cimport tools

cdef convert_1d(double *arr, int n):
	cdef Py_ssize_t i
	cdef cnp.ndarray[cnp.float64_t, ndim=1] pyarr = np.empty(n)
	for i in range(n):
		pyarr[i] = arr[i]
	return pyarr

cdef convert_2d(double **arr, int nx, int ny):
	cdef Py_ssize_t i
	cdef Py_ssize_t j
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
		cdef Py_ssize_t i_
		for i_ in range(self.argc):
			self.argv[i_] = arr[i_]

		self.Nrlk = 0
		self.Nwave = 0

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef set_wavetable(self, cnp.ndarray[double, ndim=1, mode="c"] wavetable):
		cdef Py_ssize_t i
		cdef double sigma_sq
		cdef double fact
		self.Nwave = wavetable.size
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
		cdef rh.mySpectrum spec

		spec = rh.rhf1d(self.argc, self.argv, Ndep,
				 &scale[0], &temp[0], &ne[0], &vz[0], &vmic[0],
				 &mag[0], &gamma[0], &chi[0],
				 &nH[0,0], atm_scale, &self.rlk_lines,
				 self.wavetable_ptr, self.Nwave)

		lam = convert_1d(spec.lam, spec.nlw)
		sI = convert_1d(spec.sI, spec.nlw)
		sQ, sU, sV = None, None, None
		if spec.stokes:
			sQ = convert_1d(spec.sQ, spec.nlw)
			sU = convert_1d(spec.sU, spec.nlw)
			sV = convert_1d(spec.sV, spec.nlw)
		J = convert_2d(spec.J, spec.nlw, spec.Nrays)

		return Spectrum(spec.nlw, lam, sI, sQ, sU, sV, J, None, spec.stokes)

	cpdef read_RLK_lines(self):
		self.rlk_lines = rh.get_RLK_lines(self.argc, self.argv)
		self.Nrlk = self.rlk_lines.Nrlk

	def get_Nrlk(self):
		return self.Nrlk