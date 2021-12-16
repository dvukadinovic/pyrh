# distutils: language = c
# distutils: sources = [rh/rhf1d/pyrh_compute1dray.c]
# distutils: include_dirs = rh/

cimport rh
import ctypes
import numpy as np
cimport numpy as cnp
import cython
from libc.stdlib cimport malloc
import time

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

cdef double* pyarray2double_1d(arr):
	cdef double *rh_arr
	nz = len(arr)
	rh_arr = <double *>malloc(nz*cython.sizeof(double))
	for i in range(nz):
		rh_arr[i] = arr[i]

	return rh_arr

cdef double** pyarray2double_2d(arr, Nrows, Ncols):
	cdef double **rh_matrix = rh.matrix_double(Nrows,Ncols)
	for i in range(Nrows):
		for j in range(Ncols):
			rh_matrix[i][j] = arr[i,j]

	return rh_matrix

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

# def grabnum(int[:] z):
# 	# cdef cnp.ndarray[int, ndim=1, mode="c"] arr
# 	# arr = np.ascontiguousarray(z, dtype=ctypes.c_int)
# 	a = rh._getnumber(&z[0])
# 	return a

cdef class RH:
	cdef rh.mySpectrum spec
	cdef rh.RLK_Line *rlk_lines
	cdef public int Nrlk

	def __init__(self):
		pass	

	cpdef rhf1d(self, argc, py_argv, scale, temp, ne, vz, vmic, 
		 	  mag, gamma, chi, nH, atm_scale):
		Ndep = len(temp)
		rh_scale = pyarray2double_1d(scale)
		rh_temp = pyarray2double_1d(temp)
		rh_ne = pyarray2double_1d(ne)
		rh_vz = pyarray2double_1d(vz)
		rh_vmic= pyarray2double_1d(vmic)
		rh_mag = pyarray2double_1d(mag)
		rh_gamma = pyarray2double_1d(gamma)
		rh_chi = pyarray2double_1d(chi)
		rh_nH = pyarray2double_2d(nH, 6, Ndep)

		py_list = py_argv.split(" ")
		py_string = [item.encode("utf-8") for item in py_list]
		arr = (ctypes.c_char_p * argc)(*py_string)
		cdef char *argv[10]
		for i_ in range(argc):
			argv[i_] = arr[i_]

		self.spec = rh.rhf1d(argc, argv, Ndep,
				 rh_scale, rh_temp, rh_ne, rh_vz, rh_vmic,
				 rh_mag, rh_gamma, rh_chi, 
				 rh_nH, atm_scale, 0)

		lam = convert_1d(self.spec.lam, self.spec.nlw)
		sI = convert_1d(self.spec.sI, self.spec.nlw)
		sQ, sU, sV = None, None, None
		if self.spec.stokes:
			sQ = convert_1d(self.spec.sQ, self.spec.nlw)
			sU = convert_1d(self.spec.sU, self.spec.nlw)
			sV = convert_1d(self.spec.sV, self.spec.nlw)
		J = convert_2d(self.spec.J, self.spec.nlw, self.spec.Nrays)
		
		return Spectrum(self.spec.nlw, lam, sI, sQ, sU, sV, J, None, self.spec.stokes)
	
	cpdef read_RLK_lines(self, argc, py_argv):
		py_list = py_argv.split(" ")
		py_string = [item.encode("utf-8") for item in py_list]
		arr = (ctypes.c_char_p * argc)(*py_string)
		cdef char *argv[10]
		for i_ in range(argc):
			argv[i_] = arr[i_]

		self.Nrlk = rh.get_RLK_lines(argc, argv, self.rlk_lines)

# ToDo:
#
#   -- compare the speed with and without writting to the disk
#
#   -- read input parameters
#   -- forward them to the rhf1d() and solveray() as
#      an InputData type
#   -- check if the spectrum is good
#   -- repeat this for Atoms
#   -- check consistensy
#   -- repeat this for Molecules
#   -- check consistensy
#   -- input wavelengts (not from file)

def solveray(argc, py_argv, 
	 		scale, temp, ne, 
	 		vz, vmic, 
	 		mag, gamma, chi,
			nH, atm_scale):
	Ndep = len(temp)
	rh_scale = pyarray2double_1d(scale)
	rh_temp = pyarray2double_1d(temp)
	rh_ne = pyarray2double_1d(ne)
	rh_vz = pyarray2double_1d(vz)
	rh_vmic= pyarray2double_1d(vmic)
	rh_mag = pyarray2double_1d(mag)
	rh_gamma = pyarray2double_1d(gamma)
	rh_chi = pyarray2double_1d(chi)
	rh_nH = pyarray2double_2d(nH, 6, Ndep)

	py_list = py_argv.split(" ")
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	cdef char *argv[10]
	for i_ in range(argc):
		argv[i_] = arr[i_]

	# cdef rh.mySpectrum spec
	# spec = rh._solveray(argc, argv, Ndep,
	# 				rh_scale, rh_temp, rh_ne, rh_vz, rh_vmic,
	# 		 		rh_mag, rh_gamma, rh_chi,
	# 		  		rh_nH, 1.0, atm_scale)

	# lam = convert_1d(spec.lam, spec.nlw)
	# sI = convert_1d(spec.sI, spec.nlw)
	# sQ, sU, sV = None, None, None
	# if spec.stokes:
	# 	sQ = convert_1d(spec.sQ, spec.nlw)
	# 	sU = convert_1d(spec.sU, spec.nlw)
	# 	sV = convert_1d(spec.sV, spec.nlw)
	# J = convert_2d(spec.J, spec.nlw, spec.Nrays)
	
	# return Spectrum(spec.nlw, lam, sI, sQ, sU, sV, J, None, spec.stokes)
	return 1