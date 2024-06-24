# distutils: language = c
# distutils: sources = [rh/rhf1d/pyrh_compute1dray.c]
# distutils: include_dirs = rh/

__version__ = 0.21

cimport rh
import numpy as np
cimport numpy as cnp
import cython
import ctypes
import sys

from libc.stdlib cimport malloc, free

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
	cdef cnp.ndarray[cnp.float64_t, ndim=2] pyarr = np.zeros((nx, ny))
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

# add OF table as input to rhf1d()

cdef void string2pointer(string, char* c_char_pointer):
	py_list = string.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		c_char_pointer[i_] = arr[i_]

cdef class RH:
	cdef char* cwd[160]

	# cdef int Nrlk
	# cdef rh.myRLK_Line rlk_lines

	# cdef Py_ssize_t Nwave
	# cdef array[double] wavetable
	# cdef double[::1] wavetable

	def __init__(self, cwd):
		# convert cwd to the c char pointer
		py_list = cwd.split(" ")
		argc = len(py_list)
		py_string = [item.encode("utf-8") for item in py_list]
		arr = (ctypes.c_char_p * argc)(*py_string)
		for i_ in range(argc):
			self.cwd[i_] = arr[i_]
		# string2pointer(cwd, self.cwd[0])

@cython.boundscheck(False)
@cython.wraparound(False)
def get_ne_from_nH(cwd,
					 int atm_scale,
					 cnp.ndarray[double, ndim=1, mode="c"] scale,
					 cnp.ndarray[double, ndim=1, mode="c"] temperature,
					 cnp.ndarray[double, ndim=1, mode="c"] nH,
					 cnp.ndarray[double, ndim=1, mode="c"] ne):
	cdef int Ndep = scale.size

	cdef char* argv[140]

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	rh.get_ne_from_nH(argv[0], 
			  atm_scale, 
			  Ndep, 
			  &scale[0], 
			  &temperature[0], 
			  &nH[0], 
			  &ne[0])

@cython.boundscheck(False)
@cython.wraparound(False)
def hse(cwd,
		int atm_scale,
		cnp.ndarray[double, ndim=1, mode="c"] scale,
		cnp.ndarray[double, ndim=1, mode="c"] temp,
		double pg_top,
		cnp.ndarray[double, ndim=1, mode="c"] fudge_wave=None,
		cnp.ndarray[double, ndim=2, mode="c"] fudge_value=None,
		full_output=False):
	cdef int Ndep = scale.size
	cdef int fudge_num

	cdef char* argv[140]

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	if (fudge_wave is None) or (fudge_value is None):
		myPops = rh.hse(argv[0], Ndep, pg_top,
						&scale[0], &temp[0],
						atm_scale,
						0, NULL, NULL)
	else:
		fudge_num = fudge_wave.size
		myPops = rh.hse(argv[0], Ndep, pg_top,
					 &scale[0], &temp[0],
					 atm_scale,
					 fudge_num, &fudge_wave[0], &fudge_value[0,0])

	ne = convert_1d(myPops.ne, Ndep)
	nHtot = convert_1d(myPops.nHtot, Ndep)
	# nH = convert_2d(myPops.nH, 6, Ndep)

	if full_output:
		rho = convert_1d(myPops.rho, Ndep)
		pg = convert_1d(myPops.pg, Ndep)

		return ne, nHtot, rho, pg

	return ne, nHtot

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef get_tau(cwd,
			  double mu,
			  int atm_scale,
			  cnp.ndarray[double, ndim=1, mode="c"] scale,
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

	cdef cnp.ndarray[double, ndim=1, mode="c"] tau = np.empty(Ndep)

	rh.get_tau(argv[0], mu, Ndep, &tau[0],
			 &scale[0], &atmosphere[1,0], 
			 &atmosphere[2,0], &atmosphere[3,0], &atmosphere[4,0],
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
				cnp.ndarray[double, ndim=1, mode="c"] lam_values,
				get_atomic_rfs):
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

	if get_atomic_rfs:
		rh_get_atomic_rfs = 1
	else:
		rh_get_atomic_rfs = 0

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
			 rh_get_atomic_rfs,
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

	# if (Nloggf!=0) or (Nlam!=0):
	if get_atomic_rfs:
		rf = convert_2d(spec.rfs, spec.nlw, Nloggf+Nlam)

		return sI, sQ, sU, sV, lam, rf

	# J = convert_2d(spec.J, spec.nlw, spec.Nrays)

	# Nlam = len(lam)
	return sI, sQ, sU, sV, lam

	# cpdef read_RLK_lines(self):
	# 	self.rlk_lines = rh.get_RLK_lines(self.argc, self.argv)
	# 	self.Nrlk = self.rlk_lines.Nrlk

	# def get_Nrlk(self):
	# 	return self.Nrlk

