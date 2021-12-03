# distutils: language = c
# distutils: sources = rh/rhf1d/dv_rhf1d.c
# distutils: include_dirs = rh

cimport rh
import ctypes
import numpy as np
cimport numpy as cnp
import cython
from libc.stdlib cimport malloc

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

cdef double* pyarray_to_double_1d(arr):
	cdef double *rh_arr
	nz = len(arr)
	rh_arr = <double *>malloc(nz*cython.sizeof(double))
	for i in range(nz):
		rh_arr[i] = arr[i]

	return rh_arr

cdef double** pyarray_to_double_2d(arr, Nrows, Ncols):
	cdef double **rh_matrix = rh.matrix_double(Nrows,Ncols)
	for i in range(Nrows):
		for j in range(Ncols):
			rh_matrix[i][j] = arr[i,j]

	return rh_matrix

class Spectrum(object):
	def __init__(self, nlw, lam, sI, sQ, sU, sV, J, J20, stokes):
		self.nlw = nlw
		self.lam = lam
		self.I = sI
		self.Q = sQ
		self.U = sU
		self.V = sV
		self.J = J
		self.J20 = J20
		if stokes:
			self.stokes = True
		else:
			self.stokes = False

def py_rhf1d(argc, py_argv, 
		scale, temp, ne, vz, vmic,
		mag, gamma, chi, nH, atm_scale):
	Ndep = len(temp)
	rh_scale = pyarray_to_double_1d(scale)
	rh_temp = pyarray_to_double_1d(temp)
	rh_ne = pyarray_to_double_1d(ne)
	rh_vz = pyarray_to_double_1d(vz)
	rh_vmic= pyarray_to_double_1d(vmic)
	rh_mag = pyarray_to_double_1d(mag)
	rh_gamma = pyarray_to_double_1d(gamma)
	rh_chi = pyarray_to_double_1d(chi)
	rh_nH = pyarray_to_double_2d(nH, 6, Ndep)

	py_list = py_argv.split(" ")
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	cdef char *argv[10]
	for i_ in range(argc):
		argv[i_] = arr[i_]

	cdef rh.mySpectrum spec
	spec = rh.rhf1d(argc, argv, Ndep,
				rh_scale, rh_temp, rh_ne, rh_vz, rh_vmic,
				rh_mag, rh_gamma, rh_chi, rh_nH, atm_scale)

	lam = convert_1d(spec.lam, spec.nlw)
	sI = convert_2d(spec.sI, spec.nlw, spec.Nrays)
	if spec.stokes:
		sQ = convert_2d(spec.sQ, spec.nlw, spec.Nrays)
		sU = convert_2d(spec.sU, spec.nlw, spec.Nrays)
		sV = convert_2d(spec.sV, spec.nlw, spec.Nrays)
	else:
		sQ = None
		sV = None
		sU = None
	J = convert_2d(spec.J, spec.nlw, spec.Nrays)

	return Spectrum(spec.nlw, lam, sI, sQ, sU, sV, J, None, spec.stokes)

def read_input(argc, py_argv):
	py_list = py_argv.split(" ")
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	cdef char *argv[10]
	for i_ in range(argc):
		argv[i_] = arr[i_]

	InputData = rh.readMe(argc, argv)

	return InputData