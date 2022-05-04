import numpy as np
cimport numpy as cnp
import cython
from libc.stdlib cimport malloc

cimport rh

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

"""
Maybe will be usefull sometime later
"""
cdef class ArrayWrapper:
	"""Wrap an array allocated in C"""
	cdef void* data_ptr
	cdef int nx, ny

	cdef set_data(self, int nx, void* data_ptr):
		""" Set the data of the array
		This cannot be done in the constructor as it must receive C-level
		arguments.

		Parameters:
		-----------
		nx: int
			Number of image rows
		data_ptr: void*
			Pointer to the data
		"""
		self.data_ptr = data_ptr
		self.nx = nx

	cdef as_ndarray(self, int nx, void* data_ptr):
		"""Create an `ndarray` that doesn't own the memory, we do."""
		cdef cnp.npy_intp shape[1]
		cdef cnp.ndarray ndarray

		self.set_data(nx, data_ptr)

		shape[0] = self.nx

		# Create a 2D array, of length `nx*ny/2+1`
		ndarray = cnp.PyArray_SimpleNewFromData(1, shape, cnp.NPY_FLOAT64, self.data_ptr)
		ndarray.base = <PyObject*> self

		# without this, data would be cleaned up right away
		Py_INCREF(self)
		return ndarray

	def __dealloc__(self):
		""" Frees the array. This is called by Python when all the
		references to the object are gone. """
		# print("Deallocating array")
		free(<void*>self.data_ptr)