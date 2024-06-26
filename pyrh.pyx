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
		double pg_top=0.1,
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

	cdef cnp.ndarray[cnp.float64_t, ndim=1] ne = np.empty(Ndep)
	cdef cnp.ndarray[cnp.float64_t, ndim=1] nHtot = np.empty(Ndep)
	cdef cnp.ndarray[cnp.float64_t, ndim=1] rho = np.empty(Ndep)
	cdef cnp.ndarray[cnp.float64_t, ndim=1] pg = np.empty(Ndep)
	pg[0] = pg_top

	if (fudge_wave is None) or (fudge_value is None):
		rh.hse(argv[0], Ndep,
				&scale[0], &temp[0],
				&ne[0], &nHtot[0], &rho[0], &pg[0],
				atm_scale,
				0, NULL, NULL)
	else:
		fudge_num = fudge_wave.size
		rh.hse(argv[0], Ndep,
				&scale[0], &temp[0],
				&ne[0], &nHtot[0], &rho[0], &pg[0],
				atm_scale,
				fudge_num, &fudge_wave[0], &fudge_value[0,0])

	if full_output:
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
def compute1d(cwd,
				double mu,
				int atm_scale,
				cnp.ndarray[double, ndim=2, mode="c"] atmosphere,
				cnp.ndarray[double, ndim=1, mode="c"] wave,
				cnp.ndarray[int, ndim=1, mode="c"] loggf_ids=None,
				cnp.ndarray[double, ndim=1, mode="c"] loggf_values=None,
				cnp.ndarray[int, ndim=1, mode="c"] lam_ids=None,
				cnp.ndarray[double, ndim=1, mode="c"] lam_values=None,
				cnp.ndarray[double, ndim=1, mode="c"] fudge_wave=None,
				cnp.ndarray[double, ndim=2, mode="c"] fudge_value=None,
				get_atomic_rfs=False):
	cdef int Ndep = atmosphere.shape[1]
	cdef int Nwave = wave.size
	
	#--- fudge pointers
	cdef int fudge_num = 0
	cdef double* fudge_wave_ptr = NULL
	cdef double* fudge_value_ptr = NULL
	
	if (fudge_wave is not None) or (fudge_value is not None):
		fudge_num = fudge_wave.size
		fudge_wave_ptr = &fudge_wave[0]
		fudge_value_ptr = &fudge_value[0,0]

	#--- log(gf) pointers
	cdef int Nloggf = 0
	cdef int* loggf_ids_ptr = NULL
	cdef double* loggf_values_ptr = NULL

	if (loggf_ids is not None) and (loggf_values is not None):
		Nloggf = loggf_ids.shape[0]
		if (Nloggf!=loggf_values.shape[0]):
			print("\n  pyrh: Different number of loggf_ids and loggf_values.\n")
			sys.exit()
		loggf_ids_ptr = &loggf_ids[0]
		loggf_values_ptr = &loggf_values[0]

	#--- lambda pointers
	cdef int Nlam = 0
	cdef int* lam_ids_ptr = NULL
	cdef double* lam_values_ptr = NULL

	if (lam_ids is not None) and (lam_values is not None):
		Nlam = lam_ids.shape[0]
		if (Nlam!=lam_values.shape[0]):
			print("\n  pyrh: Different number of lam_ids and lam_values.\n")
			sys.exit()
		lam_ids_ptr = &lam_ids[0]
		lam_values_ptr = &lam_values[0]

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
			 fudge_num, fudge_wave_ptr, fudge_value_ptr,
			 Nloggf, loggf_ids_ptr, loggf_values_ptr,
			 Nlam, lam_ids_ptr, lam_values_ptr,
			 rh_get_atomic_rfs,
			 0, argv[0])
			 # &self.wavetable[0], self.Nwave)

	lam = np.asarray(<cnp.float64_t[:spec.nlw]> spec.lam)
	sI = np.asarray(<cnp.float64_t[:spec.nlw]> spec.sI)

	sQ, sU, sV = None, None, None
	if spec.stokes:
		sQ = np.asarray(<cnp.float64_t[:spec.nlw]> spec.sQ)
		sU = np.asarray(<cnp.float64_t[:spec.nlw]> spec.sU)
		sV = np.asarray(<cnp.float64_t[:spec.nlw]> spec.sV)

	if get_atomic_rfs:
		rf = convert_2d(spec.rfs, spec.nlw, Nloggf+Nlam)

		return sI, sQ, sU, sV, lam, rf

	# J = convert_2d(spec.J, spec.nlw, spec.Nrays)

	return sI, sQ, sU, sV, lam